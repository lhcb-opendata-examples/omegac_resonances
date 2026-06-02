import argparse
from pathlib import Path

import matplotlib
import mplhep as hep
import numpy as np
import ROOT
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from skim import collect_files

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
hep.style.use("LHCb2")
plt.rcParams.update(
    {
        "font.size": 11,
        "axes.labelsize": 14,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 12,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "xtick.major.size": 4.5,
        "ytick.major.size": 4.5,
        "xtick.minor.width": 0.45,
        "ytick.minor.width": 0.45,
        "xtick.minor.size": 3.0,
        "ytick.minor.size": 3.0,
        "lines.linewidth": 0.9,
    }
)


def add_lhcb_label(ax):
    ax.text(
        0.04,
        0.95,
        "LHCb Open Data",
        transform=ax.transAxes,
        fontsize=16,
        fontweight="bold",
        color="#1f3b73",
        ha="left",
        va="top",
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#1f3b73", linewidth=1.0, alpha=0.96),
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build skim plots from skimmed ROOT files."
    )
    parser.add_argument(
        "skim_source",
        nargs="+",
        default=Path(__file__).resolve().parent.parent / "skimmed",
        help="Skimmed ROOT file(s) or a directory containing skimmed ROOT files.",
    )
    parser.add_argument(
        "--max-files",
        type=int,
        default=0,
        help="Limit the number of skimmed files that are read. 0 means all files.",
    )
    parser.add_argument(
        "--plot-output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "plots" / "skim",
        help="Directory where the plots will be written.",
    )
    return parser.parse_args()


def save_histogram(mass_values, output_path: Path, bins: int, mass_range, xlabel: str) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    hist, edges = np.histogram(mass_values, bins=bins, range=mass_range)
    centers = 0.5 * (edges[:-1] + edges[1:])
    errors = np.sqrt(hist)
    errors[errors == 0] = 1.0
    bin_width = edges[1] - edges[0]

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    ax.step(centers, hist, where="mid", color="black", linewidth=0.9)
    ax.errorbar(
        centers,
        hist,
        errors,
        linestyle="",
        marker=None,
        color="black",
        elinewidth=0.4,
        capsize=0,
        alpha=0.75,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(f"Candidates / {bin_width:.2f} MeV")
    ax.set_xlim(mass_range)
    ax.tick_params(axis="both", labelsize=14, width=0.65, length=4.5)
    ax.grid(True, which="both", alpha=0.2, linewidth=0.5)
    add_lhcb_label(ax)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    if len(args.skim_source) == 1:
        skim_files = collect_files(str(args.skim_source[0]), args.max_files)
    else:
        skim_files = [str(path) for path in args.skim_source]
        if args.max_files > 0:
            skim_files = skim_files[: args.max_files]
    print(f"Reading {len(skim_files)} skimmed file(s)")
    for file_name in skim_files:
        print(file_name)

    df = ROOT.RDataFrame("DecayTree/DecayTree", skim_files)
    masses = df.AsNumpy(columns=["Lambda_cplus_M", "Omega_cst0_M"])
    print(f"Entries after skim: {len(masses['Lambda_cplus_M'])}")

    print(f"Writing plots to {args.plot_output_dir}")
    save_histogram(
        np.asarray(masses["Lambda_cplus_M"], dtype=np.float64),
        args.plot_output_dir / "xic_mass.png",
        70,
        (2430, 2520),
        r"$M(\Xi_c^{+})$ [MeV]",
    )
    save_histogram(
        np.asarray(masses["Omega_cst0_M"], dtype=np.float64),
        args.plot_output_dir / "omega_c_mass.png",
        70,
        (2900, 3500),
        r"$M(\Omega_c^{0})$ [MeV]",
    )
    print("Done.")


if __name__ == "__main__":
    main()
