#!/usr/bin/env python3

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from urllib.parse import urlparse


DEFAULT_INPUT = (
    "root://eospublic.cern.ch//eos/opendata/lhcb/upload/"
    "opendata-lhcb-ntupling-service/analysis-productions/"
    "merge-requests/4413/outputs/real-production/"
)

DEFAULT_REMOTE_OUTPUT = (
    "root://eoslhcb.cern.ch//eos/opendata/lhcb/upload/"
    "example_omega_c_resonances/"
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Admin helper: run the local skim script, then upload completed "
            "skimmed ROOT files to EOS with xrdcp."
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        default=DEFAULT_INPUT,
        help="Input ROOT file, local directory, or xrootd directory passed to src/skim.py.",
    )
    parser.add_argument(
        "--skim-script",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "src" / "skim.py",
        help="Path to the normal user-facing skim.py script.",
    )
    parser.add_argument(
        "--local-output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "skimmed",
        help="Local directory where skim.py writes skimmed ROOT files before upload.",
    )
    parser.add_argument(
        "--remote-output-dir",
        default=DEFAULT_REMOTE_OUTPUT,
        help="Remote EOS/XRootD directory where skimmed ROOT files are uploaded.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Forwarded to skim.py. Use 0 for all files.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Pass -f to xrdcp and overwrite existing remote files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands but do not run skim.py or xrdcp.",
    )
    return parser.parse_args()


def run_command(command: list[str], dry_run: bool = False) -> None:
    print("+", " ".join(command))
    if dry_run:
        return
    subprocess.run(command, check=True)


def normalize_xrootd_path(path: str) -> str:
    """
    Convert root://host//eos/... into host and /eos/... for xrdfs.
    """
    parsed = urlparse(path)
    host = parsed.netloc
    eos_path = parsed.path

    while eos_path.startswith("//"):
        eos_path = eos_path[1:]

    if not host or not eos_path.startswith("/"):
        raise ValueError(f"Invalid XRootD path: {path}")

    return host, eos_path


def ensure_remote_directory(remote_dir: str, dry_run: bool = False) -> None:
    if shutil.which("xrdfs") is None:
        raise RuntimeError("xrdfs is required but was not found in PATH")

    host, eos_path = normalize_xrootd_path(remote_dir)

    command = ["xrdfs", host, "mkdir", "-p", eos_path]
    run_command(command, dry_run=dry_run)


def remote_file_url(remote_dir: str, local_file: Path) -> str:
    return f"{remote_dir.rstrip('/')}/{local_file.name}"


def upload_file(local_file: Path, remote_dir: str, overwrite: bool, dry_run: bool = False) -> None:
    if shutil.which("xrdcp") is None:
        raise RuntimeError("xrdcp is required but was not found in PATH")

    if not local_file.exists():
        raise FileNotFoundError(f"Local file does not exist: {local_file}")

    destination = remote_file_url(remote_dir, local_file)

    command = ["xrdcp"]
    if overwrite:
        command.append("-f")
    command.extend([str(local_file), destination])

    run_command(command, dry_run=dry_run)


def list_local_skim_files(local_output_dir: Path) -> list[Path]:
    return sorted(local_output_dir.glob("*.skim.root"))


def main() -> None:
    args = parse_args()

    if not args.skim_script.exists():
        raise FileNotFoundError(f"Could not find skim script: {args.skim_script}")

    before = set(list_local_skim_files(args.local_output_dir))

    skim_command = [
        sys.executable,
        str(args.skim_script),
        args.input,
        "--output-dir",
        str(args.local_output_dir),
    ]

    if args.max_files > 0:
        skim_command.extend(["--max-files", str(args.max_files)])

    print("Running local skim production...")
    run_command(skim_command, dry_run=args.dry_run)

    after = set(list_local_skim_files(args.local_output_dir))
    produced_files = sorted(after - before)

    # If files were overwritten/reproduced with the same names, `after - before`
    # would be empty. In that case, upload all skimmed files from the output dir.
    if not produced_files:
        produced_files = list_local_skim_files(args.local_output_dir)

    if not produced_files:
        raise RuntimeError(f"No skimmed files found in {args.local_output_dir}")

    print(f"Preparing remote directory: {args.remote_output_dir}")
    #ensure_remote_directory(args.remote_output_dir, dry_run=args.dry_run)

    print(f"Uploading {len(produced_files)} skimmed file(s)...")
    for local_file in produced_files:
        upload_file(
            local_file=local_file,
            remote_dir=args.remote_output_dir,
            overwrite=args.overwrite,
            dry_run=args.dry_run,
        )

    print("Done.")
    print("Remote output directory:")
    print(args.remote_output_dir)


if __name__ == "__main__":
    main()
