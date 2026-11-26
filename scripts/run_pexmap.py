"""
Run PEXMap from mass spectra in one command.

    python scripts/run_pexmap.py --spectra sample1.mzML sample2.mzML \
        --database data/octamerDB.pkl --sequences data/transcript_sequences.pkl \
        --output my_results [--engine sage|comet|msfragger]

What it does:
  1. Builds the search database (FASTA) from PEXMap's own reference isoforms, so
     every peptide the search can identify is also present in the PEXMap reference,
     plus optional contaminant sequences and decoys.
  2. Searches the spectra with the chosen search engine (default: Sage).
  3. Step 0  prepare_peptides.py   peptide-level FDR control
  4. Step 1  generate_kmers.py     peptides -> k-mers
  5. Step 2  annotate_peptides.py  gene / transcript / exon / junction mapping
  6. Step 3  verify_mapping.py     exact-substring verification

Search engines (installed separately; PEXMap does not identify spectra itself):
  sage       https://github.com/lazear/sage            open source; decoys and
             peptide-level q-values are produced by Sage and used directly.
  comet      https://uwpr.github.io/Comet/             open source; Comet adds decoys,
             PSMs are written as a Percolator (.pin) file and the FDR is estimated
             in Step 0 by target-decoy competition on Xcorr.
  msfragger  https://msfragger.nesvilab.org/           free for academic use (licence
             required); needs Java. Reversed decoys are added to the FASTA, PSMs are
             written as .pin and the FDR is estimated in Step 0 on the hyperscore.

For comet and msfragger the parameter file is generated from the engine's own
default template (comet -p / MSFragger --config), and only the settings below are
changed, so it matches the installed version. Use --engine_params to supply your
own parameter file instead.

Use --dry_run to prepare everything and print the commands without running them,
and --search_results to rerun the PEXMap steps on existing search output.
"""

import argparse
import glob
import json
import os
import pickle
import re
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(
    description="Search mass spectra and run the full PEXMap pipeline.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)
inp = parser.add_argument_group("input and output")
inp.add_argument("--spectra", nargs="+", help="one or more mzML files")
inp.add_argument("--output", required=True, help="output folder (created if needed)")
inp.add_argument("--database", required=True, help="PEXMap reference database (octamerDB.pkl)")
inp.add_argument("--sequences", required=True,
                 help="transcript sequences from extract_sequences.py (used for the search "
                      "database and for Step 3)")
inp.add_argument("--contaminants", help="optional FASTA of contaminant proteins (e.g. cRAP); "
                                        "their accessions are prefixed with CON__")
inp.add_argument("--organism", default="human", help="organism filter for Step 2 (default: human)")

eng = parser.add_argument_group("search engine")
eng.add_argument("--engine", choices=["sage", "comet", "msfragger"], default="sage",
                 help="search engine to use (default: sage)")
eng.add_argument("--engine_path",
                 help="path to the engine: the sage or comet program, or MSFragger.jar "
                      "(default: 'sage' / 'comet' on the PATH; required for msfragger)")
eng.add_argument("--java", default="java", help="Java program for MSFragger (default: 'java' on the PATH)")
eng.add_argument("--java_memory", default="16g", help="memory for MSFragger, e.g. 16g or 32g (default: 16g)")
eng.add_argument("--engine_params",
                 help="use this parameter file for the engine as it is (Sage JSON, comet.params "
                      "or fragger.params); the search settings below are then ignored, but its "
                      "database setting is replaced by the PEXMap search database")

srch = parser.add_argument_group("search settings")
srch.add_argument("--missed_cleavages", type=int, default=2, help="default: 2")
srch.add_argument("--cleave_at", default="KR", help="residues after which the enzyme cuts (default: KR = trypsin)")
srch.add_argument("--restrict", default="P", help="no cut before this residue (default: P; use '' for none)")
srch.add_argument("--min_len", type=int, default=7, help="minimum peptide length searched (default: 7)")
srch.add_argument("--max_len", type=int, default=50, help="maximum peptide length searched (default: 50)")
srch.add_argument("--static_mods", default="C:57.0215",
                  help="fixed modifications as residue:mass, comma-separated "
                       "(default: C:57.0215 = carbamidomethyl cysteine; use '' for none)")
srch.add_argument("--variable_mods", default="M:15.9949",
                  help="variable modifications as residue:mass, comma-separated "
                       "(default: M:15.9949 = oxidised methionine; use '' for none)")
srch.add_argument("--precursor_ppm", type=float, default=10.0, help="precursor tolerance ± ppm (default: 10)")
srch.add_argument("--fragment_ppm", type=float, default=20.0, help="fragment tolerance ± ppm (default: 20)")

pex = parser.add_argument_group("PEXMap settings")
pex.add_argument("--fdr", type=float, default=0.01, help="peptide-level FDR (default: 0.01)")
pex.add_argument("--min_coverage", type=float, default=0.0, help="Step 2 coverage threshold in %% (default: 0)")
pex.add_argument("--kmer_length", type=int, default=8, help="k-mer length of the database (default: 8)")

run = parser.add_argument_group("running")
run.add_argument("--dry_run", action="store_true",
                 help="prepare the database and parameters and print the commands without running them")
run.add_argument("--search_results", "--sage_results", dest="search_results",
                 help="skip the search and use existing results: results.sage.tsv for sage, "
                      "a .pin file for comet or msfragger")
args = parser.parse_args()


# --------------------------------------------------------------------------- #
# Helpers                                                                     #
# --------------------------------------------------------------------------- #

out_dir = os.path.abspath(args.output)
os.makedirs(out_dir, exist_ok=True)
log_path = os.path.join(out_dir, "run_pexmap.log")
log = open(log_path, "a")


def say(msg=""):
    print(msg)
    log.write(msg + "\n")
    log.flush()


def fail(msg):
    say("ERROR: " + msg)
    log.close()
    sys.exit(1)


def check_file(path, what):
    if not path or not os.path.isfile(path):
        fail(f"{what} not found: {path}")


def parse_mods(text, label):
    mods = {}
    for item in [t.strip() for t in text.split(",") if t.strip()]:
        if ":" not in item:
            fail(f"{label} '{item}' should look like residue:mass, e.g. C:57.0215")
        res, mass = item.split(":", 1)
        try:
            mods.setdefault(res.strip(), []).append(float(mass))
        except ValueError:
            fail(f"{label} mass '{mass}' is not a number")
    return mods


def show(cmd):
    return " ".join(f'"{c}"' if " " in c else c for c in cmd)


def run_cmd(name, cmd, cwd=None, required_output=None, always=False):
    """Run a command, log its output and stop on failure. In a dry run only print it
    (unless always=True, used for generating default parameter templates)."""
    say(f"\n=== {name} ===")
    say("$ " + show(cmd))
    if args.dry_run and not always:
        return ""
    start = time.time()
    try:
        proc = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    except OSError as e:
        fail(f"{name} could not be started: {e}")
    say(proc.stdout.rstrip())
    say(f"({time.time() - start:.0f} s)")
    if proc.returncode != 0:
        fail(f"{name} failed (exit code {proc.returncode}); see the messages above and {log_path}")
    if required_output and not os.path.exists(required_output):
        fail(f"{name} finished but did not write {required_output}")
    return proc.stdout


def find_program(name_or_path, label, url):
    exe = shutil.which(name_or_path) or (name_or_path if os.path.isfile(name_or_path) else None)
    if exe:
        exe = os.path.abspath(exe)   # engines run inside their own output folder
    if exe is None and not args.dry_run:
        fail(f"{label} was not found ('{name_or_path}'). Install it from {url} and put it on "
             f"the PATH, or give its location with --engine_path.")
    return exe or name_or_path


def set_params(path, values, essential=()):
    """Set `key = value` lines in a Comet/MSFragger parameter file, keeping comments.
    A value of None comments the line out. Keys absent from the file are appended only
    if listed in `essential` (unknown keys can make older engine versions stop).
    Returns (appended keys, skipped keys)."""
    with open(path) as f:
        lines = f.read().splitlines()
    done = set()
    for i, line in enumerate(lines):
        m = re.match(r"^\s*([A-Za-z0-9_]+)\s*=\s*([^#]*?)(\s*#.*)?$", line)
        if m and m.group(1) in values:
            key = m.group(1)
            if values[key] is None:
                lines[i] = "# " + line
            else:
                lines[i] = f"{key} = {values[key]}" + (f"  {m.group(3).strip()}" if m.group(3) else "")
            done.add(key)
    appended = [k for k in values if k not in done and k in essential and values[k] is not None]
    skipped = [k for k in values if k not in done and k not in essential and values[k] is not None]
    for k in appended:
        lines.append(f"{k} = {values[k]}")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return appended, skipped


def get_param(path, key):
    with open(path) as f:
        for line in f:
            m = re.match(r"^\s*" + re.escape(key) + r"\s*=\s*([^#]*)", line)
            if m:
                return m.group(1).strip()
    return None


def combine_pins(pins, out_path):
    """Concatenate Percolator input files (one header; drop DefaultDirection rows)."""
    header = None
    n = 0
    with open(out_path, "w") as out:
        for p in pins:
            with open(p) as f:
                first = f.readline()
                if header is None:
                    header = first
                    out.write(first)
                for line in f:
                    if line.startswith("DefaultDirection") or not line.strip():
                        continue
                    out.write(line)
                    n += 1
    return n


# --------------------------------------------------------------------------- #
# Checks                                                                      #
# --------------------------------------------------------------------------- #

say(f"PEXMap run started {time.strftime('%Y-%m-%d %H:%M:%S')}  (engine: {args.engine})")
say(f"Output folder: {out_dir}")

check_file(args.database, "PEXMap database (--database)")
check_file(args.sequences, "transcript sequences (--sequences)")
if args.contaminants:
    check_file(args.contaminants, "contaminant FASTA (--contaminants)")
if args.engine_params:
    check_file(args.engine_params, "engine parameter file (--engine_params)")
if not 0 < args.fdr < 1:
    fail("--fdr must be between 0 and 1")

static_mods = parse_mods(args.static_mods, "static modification")
variable_mods = parse_mods(args.variable_mods, "variable modification")
for res, masses in static_mods.items():
    if len(masses) > 1:
        fail(f"only one fixed modification per residue is possible ({res})")

spectra = []
if args.search_results:
    check_file(args.search_results, "search results (--search_results)")
else:
    if not args.spectra:
        fail("give the spectra with --spectra (mzML files), or existing results with --search_results")
    for s in args.spectra:
        check_file(s, "spectrum file")
        if not s.lower().endswith((".mzml", ".mzml.gz")):
            say(f"WARNING: {s} does not end in .mzML. Convert raw instrument files to mzML first "
                f"(e.g. with ProteoWizard MSConvert).")
    spectra = [os.path.abspath(s) for s in args.spectra]

engine_dir = os.path.join(out_dir, args.engine)
os.makedirs(engine_dir, exist_ok=True)


# --------------------------------------------------------------------------- #
# 1. Search database from the PEXMap reference                                #
# --------------------------------------------------------------------------- #

fasta_path = os.path.join(out_dir, "search_database.fasta")


def write_fasta_record(out, acc, seq):
    out.write(f">{acc}\n")
    for i in range(0, len(seq), 60):
        out.write(seq[i:i + 60] + "\n")


if not args.search_results:
    say("\n=== Building the search database from the PEXMap reference ===")
    with open(args.sequences, "rb") as f:
        sequences = pickle.load(f)
    if not isinstance(sequences, dict) or not sequences:
        fail("--sequences should be the transcript sequence file written by extract_sequences.py")
    records = []
    for tid in sorted(sequences):
        seq = str(sequences[tid]).replace("*", "").strip().upper()
        if seq:
            records.append((str(tid), seq))
    n_target = len(records)
    if args.contaminants:
        acc, buf = None, []
        with open(args.contaminants) as c:
            for line in c:
                if line.startswith(">"):
                    if acc:
                        records.append((acc, "".join(buf)))
                    head = line[1:].strip()
                    acc = head.split()[0] if head else f"unnamed{len(records)}"
                    if not acc.startswith(("CON__", "contam_")):
                        acc = "CON__" + acc
                    buf = []
                elif line.strip():
                    buf.append(line.strip().upper())
            if acc:
                records.append((acc, "".join(buf)))
    n_cont = len(records) - n_target
    # Sage and Comet generate decoys themselves; MSFragger needs them in the FASTA.
    add_decoys = args.engine == "msfragger"
    with open(fasta_path, "w") as out:
        for acc, seq in records:
            write_fasta_record(out, acc, seq)
        if add_decoys:
            for acc, seq in records:
                write_fasta_record(out, "rev_" + acc, seq[::-1])
    say(f"Reference isoform sequences: {n_target}")
    say(f"Contaminant sequences: {n_cont}" + ("" if args.contaminants else
        " (none given; consider --contaminants with a cRAP FASTA)"))
    say("Decoys: " + ("reversed sequences added to the FASTA (prefix rev_)" if add_decoys else
                      f"generated by {args.engine} (reversed sequences)"))
    say(f"Written: {fasta_path}")


# --------------------------------------------------------------------------- #
# 2. Search                                                                   #
# --------------------------------------------------------------------------- #

# What Step 0 needs from each engine
if args.engine == "sage":
    step0_format = ["--format", "sage"]
    results = args.search_results or os.path.join(engine_dir, "results.sage.tsv")
else:
    score = "Xcorr" if args.engine == "comet" else "hyperscore"
    step0_format = ["--format", "pin", "--score_col", score]
    results = args.search_results or os.path.join(engine_dir, "combined.pin")

if not args.search_results:

    # ---------------- Sage ----------------
    if args.engine == "sage":
        exe = find_program(args.engine_path or "sage", "The Sage program", "https://github.com/lazear/sage/releases")
        config_path = os.path.join(engine_dir, "sage_config.json")
        if args.engine_params:
            with open(args.engine_params) as f:
                config = json.load(f)
            config.setdefault("database", {})["fasta"] = fasta_path
            config["database"]["generate_decoys"] = True
            config["output_directory"] = engine_dir
            config["mzml_paths"] = spectra
        else:
            config = {
                "database": {
                    "bucket_size": 32768,
                    "enzyme": {
                        "missed_cleavages": args.missed_cleavages,
                        "min_len": args.min_len,
                        "max_len": args.max_len,
                        "cleave_at": args.cleave_at,
                        "restrict": args.restrict if args.restrict else None,
                    },
                    "max_variable_mods": 2,
                    "decoy_tag": "rev_",
                    "generate_decoys": True,
                    "fasta": fasta_path,
                },
                "precursor_tol": {"ppm": [-args.precursor_ppm, args.precursor_ppm]},
                "fragment_tol": {"ppm": [-args.fragment_ppm, args.fragment_ppm]},
                "report_psms": 1,
                "output_directory": engine_dir,
                "mzml_paths": spectra,
            }
            if static_mods:
                config["database"]["static_mods"] = {r: m[0] for r, m in static_mods.items()}
            if variable_mods:
                config["database"]["variable_mods"] = variable_mods
        with open(config_path, "w") as f:
            json.dump(config, f, indent=2)
        say(f"\nSage configuration written: {config_path}")
        run_cmd("Search: Sage", [exe, config_path, "-o", engine_dir] + spectra, required_output=results)

    # ---------------- Comet ----------------
    elif args.engine == "comet":
        exe = find_program(args.engine_path or "comet", "The Comet program", "https://uwpr.github.io/Comet/")
        params = os.path.join(engine_dir, "comet.params")
        if args.engine_params:
            shutil.copy(args.engine_params, params)
        elif args.dry_run and not shutil.which(exe) and not os.path.isfile(exe):
            say("\n(dry run: Comet not found, so its default parameter template cannot be generated)")
            params = None
        else:
            run_cmd("Comet default parameters", [exe, "-p"], cwd=engine_dir, always=True)
            template = os.path.join(engine_dir, "comet.params.new")
            if not os.path.isfile(template):
                fail("Comet did not write comet.params.new when run with -p")
            shutil.move(template, params)
        if params:
            values = {
                "database_name": fasta_path,
                "decoy_search": "1",               # concatenated target-decoy search
                "decoy_prefix": "rev_",
                "output_percolatorfile": "1",
                "output_txtfile": "0",
                "output_pepxmlfile": "0",
                "output_sqtfile": "0",
                "output_mzidentmlfile": "0",
                "num_output_lines": "1",
            }
            if not args.engine_params:
                values.update({
                    "peptide_mass_units": "2",      # ppm
                    "allowed_missed_cleavage": str(args.missed_cleavages),
                    "fragment_bin_tol": "0.02",     # high-resolution fragments
                    "fragment_bin_offset": "0.0",
                    "theoretical_fragment_ions": "0",
                })
                if args.cleave_at != "KR" or args.restrict != "P":
                    say("WARNING: Comet selects enzymes by number from its own list; the generated "
                        "parameters keep Comet's default enzyme. For other enzymes edit comet.params "
                        "and pass it with --engine_params.")
                # tolerance keys differ between Comet versions; set whichever the template has
                if get_param(params, "peptide_mass_tolerance_upper") is not None:
                    values["peptide_mass_tolerance_upper"] = str(args.precursor_ppm)
                    values["peptide_mass_tolerance_lower"] = str(-args.precursor_ppm)
                else:
                    values["peptide_mass_tolerance"] = str(args.precursor_ppm)
                for res, mass in static_mods.items():
                    key = {"C": "add_C_cysteine", "K": "add_K_lysine"}.get(res)
                    if key:
                        values[key] = str(mass[0])
                    else:
                        say(f"WARNING: fixed modification on {res} not set for Comet; edit comet.params "
                            f"and pass it with --engine_params.")
                if "C" not in static_mods:
                    values["add_C_cysteine"] = "0.0"
                # variable mods: keep each template line's format, change only mass and residues
                vm = [(m, r) for r, masses in variable_mods.items() for m in masses]
                for i in range(1, 16):
                    key = f"variable_mod{i:02d}"
                    line = get_param(params, key)
                    if line is None:
                        continue
                    fields = line.split()
                    if i <= len(vm):
                        fields[0], fields[1] = str(vm[i - 1][0]), vm[i - 1][1]
                    else:
                        fields[0], fields[1] = "0.0", "X"
                    values[key] = " ".join(fields)
            appended, skipped = set_params(params, values,
                                           essential=("database_name", "decoy_search", "decoy_prefix",
                                                      "output_percolatorfile"))
            if appended:
                say("Note: added to comet.params (not in this Comet version's template): " + ", ".join(appended))
            if skipped:
                say("Note: not in this Comet version's parameters, left unchanged: " + ", ".join(skipped))
            say(f"Comet parameters written: {params}")
        run_cmd("Search: Comet", [exe, f"-P{params or 'comet.params'}"] + spectra, cwd=engine_dir)

    # ---------------- MSFragger ----------------
    else:
        if not args.engine_path:
            fail("for MSFragger, give the location of MSFragger.jar with --engine_path "
                 "(download it from https://msfragger.nesvilab.org/; licence required)")
        jar = args.engine_path
        if not args.dry_run:
            check_file(jar, "MSFragger.jar (--engine_path)")
        java = find_program(args.java, "Java (needed for MSFragger)", "https://adoptium.net/")
        params = os.path.join(engine_dir, "fragger.params")
        if args.engine_params:
            shutil.copy(args.engine_params, params)
        elif args.dry_run and not os.path.isfile(jar):
            say("\n(dry run: MSFragger.jar not found, so its default parameter template cannot be generated)")
            params = None
        else:
            run_cmd("MSFragger default parameters", [java, "-jar", os.path.abspath(jar), "--config"],
                    cwd=engine_dir, always=True)
            template = os.path.join(engine_dir, "closed_fragger.params")
            if not os.path.isfile(template):
                fail("MSFragger did not write closed_fragger.params when run with --config")
            shutil.move(template, params)
        if params:
            values = {
                "database_name": fasta_path,
                "output_format": "tsv_pin",
                "output_report_topN": "1",
            }
            if not args.engine_params:
                values.update({
                    "precursor_mass_lower": str(-args.precursor_ppm),
                    "precursor_mass_upper": str(args.precursor_ppm),
                    "precursor_mass_units": "1",    # ppm
                    "precursor_true_tolerance": str(args.precursor_ppm),
                    "precursor_true_units": "1",
                    "fragment_mass_tolerance": str(args.fragment_ppm),
                    "fragment_mass_units": "1",     # ppm
                    "search_enzyme_cut_1": args.cleave_at,
                    "search_enzyme_nocut_1": args.restrict if args.restrict else "",
                    "allowed_missed_cleavage_1": str(args.missed_cleavages),
                    "digest_min_length": str(args.min_len),
                    "digest_max_length": str(args.max_len),
                })
                if get_param(params, "allowed_missed_cleavage_1") is None:      # older versions
                    values["allowed_missed_cleavage"] = str(args.missed_cleavages)
                for res, mass in static_mods.items():
                    names = {"C": "add_C_cysteine", "K": "add_K_lysine", "M": "add_M_methionine"}
                    if res in names:
                        values[names[res]] = str(mass[0])
                    else:
                        say(f"WARNING: fixed modification on {res} not set for MSFragger; edit "
                            f"fragger.params and pass it with --engine_params.")
                if "C" not in static_mods:
                    values["add_C_cysteine"] = "0.0"
                vm = [(m, r) for r, masses in variable_mods.items() for m in masses]
                for i in range(1, 17):
                    key = f"variable_mod_{i:02d}"
                    existing = get_param(params, key)
                    if i <= len(vm):
                        values[key] = f"{vm[i - 1][0]} {vm[i - 1][1]} 3"
                    elif existing is not None:
                        values[key] = None        # comment out template defaults beyond ours
            appended, skipped = set_params(params, values,
                                           essential=("database_name", "output_format") +
                                           tuple(f"variable_mod_{i:02d}" for i in range(1, 17)))
            if appended:
                say("Note: added to fragger.params (not in this MSFragger version's template): " + ", ".join(appended))
            if skipped:
                say("Note: not in this MSFragger version's parameters, left unchanged: " + ", ".join(skipped))
            say(f"MSFragger parameters written: {params}")
        run_cmd("Search: MSFragger",
                [java, f"-Xmx{args.java_memory}", "-jar", os.path.abspath(jar), params or "fragger.params"] + spectra,
                cwd=engine_dir)

    # Comet and MSFragger write <name>.pin next to each mzML file: collect and combine them
    if args.engine in ("comet", "msfragger") and not args.dry_run:
        pins = []
        for s in spectra:
            base = re.sub(r"\.mzml(\.gz)?$", "", s, flags=re.I)
            cands = [base + ".pin"] + glob.glob(base + "*.pin") + \
                    [os.path.join(engine_dir, os.path.basename(base) + ".pin")]
            found = next((c for c in cands if os.path.isfile(c)), None)
            if not found:
                fail(f"no .pin output found for {s} (looked for {base}.pin)")
            dest = os.path.join(engine_dir, os.path.basename(found))
            if os.path.abspath(found) != os.path.abspath(dest):
                shutil.move(found, dest)
            pins.append(dest)
        n = combine_pins(pins, results)
        say(f"Collected {len(pins)} .pin file(s), {n} PSMs: {results}")


# --------------------------------------------------------------------------- #
# 3-6. PEXMap steps                                                           #
# --------------------------------------------------------------------------- #

py = sys.executable
peptides = os.path.join(out_dir, "peptides.txt")
kmers = os.path.join(out_dir, "kmers.txt")
annotations = os.path.join(out_dir, "annotations.tsv")
verified = os.path.join(out_dir, "annotations_verified.tsv")

s0 = run_cmd("Step 0: peptide-level FDR control",
             [py, os.path.join(HERE, "prepare_peptides.py"), "--input", results, "--output", peptides,
              "--fdr", str(args.fdr)] + step0_format, required_output=peptides)
s1 = run_cmd("Step 1: k-mers",
             [py, os.path.join(HERE, "generate_kmers.py"), "--input", peptides,
              "--output", kmers, "--kmer_length", str(args.kmer_length)], required_output=kmers)
s2 = run_cmd("Step 2: mapping to genes, transcripts, exons and junctions",
             [py, os.path.join(HERE, "annotate_peptides.py"), "--kmers", kmers,
              "--database", args.database, "--organism", args.organism,
              "--peptides", peptides, "--kmer_length", str(args.kmer_length),
              "--min_coverage", str(args.min_coverage), "--output", annotations], required_output=annotations)
s3 = run_cmd("Step 3: exact-substring verification",
             [py, os.path.join(HERE, "verify_mapping.py"), "--annotations", annotations,
              "--sequences", args.sequences, "--output", verified, "--drop_failed"], required_output=verified)


# --------------------------------------------------------------------------- #
# Final summary                                                               #
# --------------------------------------------------------------------------- #

if args.dry_run:
    say("\nDry run: the search and PEXMap steps were not run. Remove --dry_run to run them.")
    log.close()
    sys.exit(0)


def pick(text, starts):
    for line in text.splitlines():
        if line.strip().startswith(starts):
            return line.strip()
    return "(none)"


summary = [
    "PEXMap run summary",
    "==================",
    f"Search engine: {args.engine}",
    f"Spectra: {', '.join(os.path.basename(s) for s in spectra) if spectra else '(existing search results)'}",
    f"Search results: {results}",
    "",
    "Peptide identification and FDR (Step 0)",
    "  " + pick(s0, "FDR control"),
    "  " + pick(s0, "PSMs read"),
    "  " + pick(s0, "Distinct target peptides"),
    "  " + pick(s0, "Kept"),
    "",
    "Mapping (Steps 1-2)",
    "  " + pick(s1, "Processed unique peptides"),
    "  " + pick(s2, "Annotated"),
    "  " + pick(s2, "Unmapped"),
    "  " + pick(s2, "Peptides mapped to exons"),
    "  " + pick(s2, "Peptides mapped to junctions"),
    "  " + pick(s2, "Peptides with tied"),
    "",
    "Verification (Step 3)",
    "  " + (pick(s3, "Exact-substring verified") if "Exact-substring verified" in s3
            else "No annotated peptides to verify"),
    ("  " + pick(s3, "Failed verification")) if "Failed verification" in s3 else "",
    "",
    "Main result: " + verified,
    "All messages: " + log_path,
]
text = "\n".join(summary)
say("\n" + text)
with open(os.path.join(out_dir, "summary.txt"), "w") as f:
    f.write(text + "\n")
log.close()
