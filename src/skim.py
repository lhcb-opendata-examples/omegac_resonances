import argparse
import shutil
import subprocess
from pathlib import Path
from urllib.parse import urlparse

import ROOT

ROOT.ROOT.EnableImplicitMT()
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

DEFAULT_INPUT = (
    "root://eospublic.cern.ch//eos/opendata/lhcb/upload/"
    "opendata-lhcb-ntupling-service/analysis-productions/"
    "merge-requests/4413/outputs/real-production/"
)

# The DTF-based observables are intentionally left out for now.
# The branch name `Lambda_cplus_*` comes from the stripping line naming in the
# prompt-charm selection, which reconstructs both Lambda_c+ and Xi_c+ candidates.
# Here we use `Lambda_cplus_M` as the Xi_c+ mass selection branch.
# This skim is the single offline selection for the analysis, so the tighter
# Xi_c+ quality cuts that used to live in the sweight step are included here too.

# Please play around with different sets of cuts
SKIM_CUTS = " && ".join(
    [
        "(Omega_cst0_OWNPV_CHI2 / Omega_cst0_OWNPV_NDOF) < 3",
        "Lambda_cplus_M > 2400 && Lambda_cplus_M < 2520",
        "Omega_cst0_PT > 4500",
        "Kminus_0_PT > 400",
        "(Omega_cst0_ENDVERTEX_CHI2 / Omega_cst0_ENDVERTEX_NDOF) < 3",
        "(Lambda_cplus_ENDVERTEX_CHI2 / Lambda_cplus_ENDVERTEX_NDOF) < 3",
        "(Lambda_cplus_OWNPV_CHI2 / Lambda_cplus_OWNPV_NDOF) < 3",
        "(Lambda_cplus_ORIVX_CHI2 / Lambda_cplus_ORIVX_NDOF) < 3",
        "(Kminus_OWNPV_CHI2 / Kminus_OWNPV_NDOF) < 3",
        "(pplus_OWNPV_CHI2 / pplus_OWNPV_NDOF) < 3",
        "(piplus_OWNPV_CHI2 / piplus_OWNPV_NDOF) < 3",
        "(Kminus_0_OWNPV_CHI2 / Kminus_0_OWNPV_NDOF) < 3",
        "pplus_ProbNNp > 0.9",
        "Kminus_ProbNNk > 0.8",
        "piplus_ProbNNpi > 0.8",
        "(Lambda_cplus_ENDVERTEX_CHI2 / Lambda_cplus_ENDVERTEX_NDOF) < 2",
        "pplus_ProbNNghost < 0.1",
        "pplus_ProbNNp > 0.5",
        "piplus_ProbNNpi > 0.5",
        "Kminus_ProbNNk > 0.5",
        "Kminus_0_ProbNNk > 0.5",
    ]
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Read ROOT files and build the skimmed dataframe."
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=DEFAULT_INPUT,
        help="Local directory or xrootd directory containing ROOT files",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Limit the number of files to read from the directory. Use 0 for all matching files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "skimmed",
        help="Directory where the skimmed ROOT files will be written locally.",
    )
    return parser.parse_args()


def collect_files(source: str, max_files: int) -> list[str]:
    if source.startswith("root://") and source.endswith(".root") and not source.endswith("/"):
        files = [source]
    elif source.startswith("root://"):
        files = collect_xrootd_files(source)
    else:
        path = Path(source)
        if path.is_file():
            files = [str(path)]
        else:
            if not path.is_dir():
                raise NotADirectoryError(f"Input path is not a directory or ROOT file: {source}")

            files = sorted(str(p) for p in path.glob("*.root"))

    if max_files > 0:
        files = files[:max_files]

    if not files:
        raise RuntimeError(f"No ROOT files matched input: {source}")

    return files


def collect_xrootd_files(source: str) -> list[str]:
    parsed = urlparse(source)
    host = parsed.netloc
    directory = parsed.path

    if not host or not directory:
        raise ValueError(f"Invalid xrootd directory: {source}")

    if directory.startswith("//"):
        directory = directory[1:]
    elif not directory.startswith("/"):
        directory = f"/{directory}"

    if shutil.which("xrdfs") is None:
        raise RuntimeError("xrdfs is required to list files from xrootd sources")

    command = ["xrdfs", host, "ls", directory]
    result = subprocess.run(command, check=True, capture_output=True, text=True)

    files = []
    for line in result.stdout.splitlines():
        line = line.strip()
        if line.endswith(".root"):
            files.append(f"root://{host}//{line.lstrip('/')}")

    return sorted(files)


def skimmed_output_path(output_dir: Path, input_file: str) -> Path:
    name = Path(urlparse(input_file).path).name
    return output_dir / f"{Path(name).stem}.skim.root"


def skim_file(file_name: str, output_file: Path) -> int:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df = ROOT.RDataFrame("DecayTree/DecayTree", file_name)
    skimmed_df = df.Filter(SKIM_CUTS)
    count = skimmed_df.Count()
    snapshot = skimmed_df.Snapshot("DecayTree/DecayTree", str(output_file))
    ROOT.RDF.RunGraphs([count, snapshot])
    entries = count.GetValue()
    print(f"Entries after skim: {entries}")
    print(f"Wrote {output_file}")
    return entries


def main():
    args = parse_args()

    files = collect_files(args.input, args.max_files)
    print(f"Reading {len(files)} file(s)")
    for file_name in files:
        print(file_name)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Skimmed dataframe ready with cut string:\n{SKIM_CUTS}")

    for index, file_name in enumerate(files, start=1):
        print(f"Running event loop for file {index}/{len(files)}...")
        output_file = skimmed_output_path(args.output_dir, file_name)
        skim_file(file_name, output_file)
        print(f"Finished file {index}/{len(files)}")

    print("Done.")


if __name__ == "__main__":
    main()
