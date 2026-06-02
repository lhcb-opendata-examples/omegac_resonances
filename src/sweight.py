import argparse
import json
from array import array
from pathlib import Path

import numpy as np
import ROOT
import matplotlib
import zfit

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
from sweights import SWeight

from skim import collect_files


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
hep.style.use("LHCb2")
plt.rcParams.update(
    {
        "font.size": 11,
        "axes.labelsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
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

XIC_MEAN_PDG = 2467.79
XIC_FIT_LOW = 2430.0
XIC_FIT_HIGH = 2510.0
XIC_PLOT_BINS = 70


def add_lhcb_label(ax):
    ax.text(
        0.06,
        0.80,
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
        description="Fit the Xi_c mass and compute signal sweights."
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
        "--output-root",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "weighted" / "omegac2xicK.root",
        help="Path to the weighted ROOT file to write.",
    )
    parser.add_argument(
        "--plot-output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "plots" / "sweight",
        help="Directory where the diagnostic plots will be written.",
    )
    return parser.parse_args()


def read_mass_arrays(skim_source: str, max_files: int) -> dict[str, np.ndarray]:
    if len(skim_source) == 1:
        skim_files = collect_files(str(skim_source[0]), max_files)
    else:
        skim_files = [str(path) for path in skim_source]
        if max_files > 0:
            skim_files = skim_files[:max_files]
    print(f"Reading {len(skim_files)} skimmed file(s)")
    for file_name in skim_files:
        print(file_name)
    df = ROOT.RDataFrame("DecayTree/DecayTree", skim_files)
    return df.AsNumpy(columns=["Lambda_cplus_M", "Omega_cst0_M"])


def build_signal_pdf(obs, suffix: str):
    mean = zfit.Parameter(f"xic_mean_{suffix}", XIC_MEAN_PDG, 2465.0, 2471.0)
    sigma = zfit.Parameter(f"xic_sigma_{suffix}", 4.8, 1.0, 20.0)
    alpha_left = zfit.Parameter(f"xic_alpha_left_{suffix}", 2.5, 0.5, 8.0)
    n_left = zfit.Parameter(f"xic_n_left_{suffix}", 8.0, 0.5, 80.0)
    alpha_right = zfit.Parameter(f"xic_alpha_right_{suffix}", 2.5, 0.5, 8.0)
    n_right = zfit.Parameter(f"xic_n_right_{suffix}", 8.0, 0.5, 80.0)

    signal = zfit.pdf.DoubleCB(
        obs=obs,
        mu=mean,
        sigma=sigma,
        alphal=alpha_left,
        nl=n_left,
        alphar=alpha_right,
        nr=n_right,
    )
    params = {
        "mean": mean,
        "sigma": sigma,
        "alpha_left": alpha_left,
        "n_left": n_left,
        "alpha_right": alpha_right,
        "n_right": n_right,
    }
    return signal, params


def build_zfit_model(n_events: int):
    obs = zfit.Space("Lambda_cplus_M", limits=(XIC_FIT_LOW, XIC_FIT_HIGH))
    signal, signal_params = build_signal_pdf(obs, "data")
    nsig = zfit.Parameter("nsig", 0.7 * n_events, 1.0, 2.0 * n_events)
    nbkg = zfit.Parameter("nbkg", 0.3 * n_events, 1.0, 2.0 * n_events)
    slope = zfit.Parameter("lambda", -0.0005, -0.01, 0.0)

    signal_ext = signal.create_extended(nsig)
    background = zfit.pdf.Exponential(obs=obs, lam=slope)
    background_ext = background.create_extended(nbkg)
    model = zfit.pdf.SumPDF([signal_ext, background_ext])
    params = {
        "nsig": nsig,
        "nbkg": nbkg,
        "slope": slope,
    }
    params.update(signal_params)
    return obs, model, signal_ext, background_ext, params


def fit_mass(m_xi: np.ndarray):
    obs, total_model, signal_ext, background_ext, params = build_zfit_model(len(m_xi))
    data = zfit.Data.from_numpy(obs=obs, array=np.asarray(m_xi, dtype=float))
    loss = zfit.loss.ExtendedUnbinnedNLL(model=total_model, data=data)
    minimizer = zfit.minimize.Minuit()
    result = minimizer.minimize(loss)
    try:
        fit_errors = result.hesse()
    except Exception:
        fit_errors = {}
    return data, obs, total_model, signal_ext, background_ext, params, result, fit_errors


def parameter_value(param) -> float:
    return float(param.value().numpy())


def parameter_error(fit_errors, param) -> float:
    for key in (param, getattr(param, "name", None)):
        if key is None:
            continue
        try:
            return float(fit_errors[key]["error"])
        except Exception:
            pass
    return float("nan")


def save_xic_fit_plot(
    m_xi: np.ndarray,
    obs,
    signal_ext,
    background_ext,
    params,
    fit_errors,
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    y, edges = np.histogram(m_xi, bins=XIC_PLOT_BINS, density=False, range=(XIC_FIT_LOW, XIC_FIT_HIGH))
    bin_delta = edges[1] - edges[0]
    bins = 0.5 * (edges[:-1] + edges[1:])

    xs = np.linspace(XIC_FIT_LOW, XIC_FIT_HIGH, 300)
    stacked_xs = xs.reshape((-1, 1))
    stacked_bins = bins.reshape((-1, 1))
    signal = parameter_value(params["nsig"]) * np.asarray(
        signal_ext.pdf(zfit.Data.from_numpy(obs=obs, array=stacked_xs, guarantee_limits=True)),
        dtype=np.float64,
    ) * bin_delta
    background = parameter_value(params["nbkg"]) * np.asarray(
        background_ext.pdf(zfit.Data.from_numpy(obs=obs, array=stacked_xs, guarantee_limits=True)),
        dtype=np.float64,
    ) * bin_delta

    fig, (ax_main, ax_pull) = plt.subplots(
        2,
        1,
        figsize=(8.2, 5.5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0},
    )
    ax_main.errorbar(
        bins,
        y,
        np.sqrt(y),
        linestyle="",
        marker="x",
        ms=3.8,
        color="black",
        elinewidth=0.9,
        capsize=0,
        alpha=0.8,
    )
    ax_main.plot(xs, signal, label="signal", linewidth=0.9)
    ax_main.plot(xs, background, label="background", linewidth=0.9)
    ax_main.plot(xs, signal + background, label="combined", linewidth=0.9)
    ax_main.set_ylabel(f"Candidates / {bin_delta:.2f} MeV")
    ax_main.legend(fontsize=14)
    ax_main.tick_params(axis="both", labelsize=14, width=0.65, length=4.5)
    ax_main.grid(True, which="both", alpha=0.2, linewidth=0.5)
    add_lhcb_label(ax_main)

    nsig = parameter_value(params["nsig"])
    nbkg = parameter_value(params["nbkg"])
    nsig_err = parameter_error(fit_errors, params["nsig"])
    nbkg_err = parameter_error(fit_errors, params["nbkg"])
    ax_main.text(
        0.78,
        0.63,
        (
            rf"$N_\mathrm{{sig}} = {nsig:.0f} \pm {nsig_err:.0f}$"
            "\n"
            rf"$N_\mathrm{{bkg}} = {nbkg:.0f} \pm {nbkg_err:.0f}$"
        ),
        transform=ax_main.transAxes,
        ha="center",
        va="top",
        fontsize=12,
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="#1f3b73",
            linewidth=1.0,
            alpha=0.95,
        ),
    )

    model_y = (
        parameter_value(params["nsig"])
        * np.asarray(signal_ext.pdf(zfit.Data.from_numpy(obs=obs, array=stacked_bins, guarantee_limits=True)), dtype=np.float64)
        * bin_delta
        + parameter_value(params["nbkg"])
        * np.asarray(background_ext.pdf(zfit.Data.from_numpy(obs=obs, array=stacked_bins, guarantee_limits=True)), dtype=np.float64)
        * bin_delta
    )
    pull = (y - model_y) / np.where(np.sqrt(y) > 0, np.sqrt(y), 1.0)
    ax_pull.bar(
        bins,
        pull,
        width=bin_delta,
        align="center",
        color="black",
        edgecolor="black",
        linewidth=0.25,
    )
    ax_pull.axhline(0, color="gray", linestyle="-", linewidth=0.7, alpha=0.45)
    ax_pull.axhline(-3, color="red", linestyle="--", alpha=0.5)
    ax_pull.axhline(3, color="red", linestyle="--", alpha=0.5)
    ax_pull.set_ylabel("Pull [$\\sigma$]")
    ax_pull.tick_params(axis="both", labelsize=14, width=0.65, length=4.5)
    ax_pull.set_xlabel(r"$m(\Xi_c^+)$ / MeV")
    ax_pull.set_ylim((-5, 5))
    ax_pull.grid(True, which="both", alpha=0.2, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_weighted_omega_plot(m_omega: np.ndarray, weights: np.ndarray, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_range = (2970, 3500)
    bins = 70

    fig, (ax_sig, ax_bkg) = plt.subplots(1, 2, figsize=(12.5, 4.8))
    ax_sig.hist(m_omega, bins=bins, range=plot_range, weights=weights, histtype="step", color="black", label="signal", linewidth=0.9)
    ax_sig.legend(fontsize=12)
    ax_sig.tick_params(axis="both", labelsize=14, width=0.65, length=4.5)
    ax_sig.set_xlabel(r"$m(\Omega_c^0)$")
    ax_sig.set_ylabel(f"Weighted candidates / {(plot_range[1] - plot_range[0]) / bins:.2f} MeV")
    ax_sig.set_xlim(plot_range)
    ax_sig.grid(True, which="both", alpha=0.2, linewidth=0.5)
    add_lhcb_label(ax_sig)
    ax_bkg.hist(m_omega, bins=bins, range=plot_range, weights=1 - weights, histtype="step", color="tab:red", label="background", linewidth=0.9)
    ax_bkg.legend(fontsize=12)
    ax_bkg.tick_params(axis="both", labelsize=14, width=0.65, length=4.5)
    ax_bkg.set_xlabel(r"$m(\Omega_c^0)$")
    ax_bkg.set_ylabel(f"Weighted candidates / {(plot_range[1] - plot_range[0]) / bins:.2f} MeV")
    ax_bkg.set_xlim(plot_range)
    ax_bkg.grid(True, which="both", alpha=0.2, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(output_dir / "omega_weighted.png", dpi=200)
    plt.close(fig)


def write_weighted_root(output_file: Path, masses: dict[str, np.ndarray], weights: np.ndarray) -> None:
    output_file.parent.mkdir(parents=True, exist_ok=True)
    m_xi = np.asarray(masses["Lambda_cplus_M"], dtype=np.float64)
    m_omega = np.asarray(masses["Omega_cst0_M"], dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    background = np.asarray(1.0 - weights, dtype=np.float64)

    output = ROOT.TFile(str(output_file), "RECREATE")
    directory = output.mkdir("DecayTree")
    directory.cd()
    tree = ROOT.TTree("DecayTree", "DecayTree")

    branch_m_xi = array("d", [0.0])
    branch_m_omega = array("d", [0.0])
    branch_weight = array("d", [0.0])
    branch_background = array("d", [0.0])

    tree.Branch("Lambda_cplus_M", branch_m_xi, "Lambda_cplus_M/D")
    tree.Branch("Omega_cst0_M", branch_m_omega, "Omega_cst0_M/D")
    tree.Branch("sweight_signal", branch_weight, "sweight_signal/D")
    tree.Branch("sweight_background", branch_background, "sweight_background/D")

    for xi_mass, omega_mass, sw_sig, sw_bkg in zip(m_xi, m_omega, weights, background):
        branch_m_xi[0] = float(xi_mass)
        branch_m_omega[0] = float(omega_mass)
        branch_weight[0] = float(sw_sig)
        branch_background[0] = float(sw_bkg)
        tree.Fill()

    tree.Write()
    output.Close()


def main():
    args = parse_args()
    output_root = args.output_root or (Path(__file__).resolve().parent.parent / "weighted" / "omegac2xicK.root")

    masses = read_mass_arrays(args.skim_source, args.max_files)
    m_xi = np.asarray(masses["Lambda_cplus_M"], dtype=np.float64)
    m_omega = np.asarray(masses["Omega_cst0_M"], dtype=np.float64)

    print(f"Fitting Xi_c mass with {len(m_xi)} entries")
    data, obs, total_model, signal_ext, background_ext, params, fit_result, fit_errors = fit_mass(m_xi)
    syield = parameter_value(params["nsig"])
    byield = parameter_value(params["nbkg"])

    print("Building sweights")
    def spdf(m):
        values = np.asarray(m, dtype=np.float64)
        data = zfit.Data.from_numpy(obs=obs, array=values.reshape((-1, 1)), guarantee_limits=True)
        return np.asarray(signal_ext.pdf(data), dtype=np.float64)

    def bpdf(m):
        values = np.asarray(m, dtype=np.float64)
        data = zfit.Data.from_numpy(obs=obs, array=values.reshape((-1, 1)), guarantee_limits=True)
        return np.asarray(background_ext.pdf(data), dtype=np.float64)

    sweight = SWeight(
        data=m_xi,
        pdfs=[spdf, bpdf],
        yields=[syield, byield],
        discvarranges=[(2430.0, 2510.0)],
        verbose=True,
        checks=True,
    )
    weights = sweight(m_xi)

    print(f"Writing diagnostic plots to {args.plot_output_dir}")
    args.plot_output_dir.mkdir(parents=True, exist_ok=True)
    save_xic_fit_plot(
        m_xi,
        obs,
        signal_ext,
        background_ext,
        params,
        fit_errors,
        args.plot_output_dir / "xic_fit.png",
    )
    save_weighted_omega_plot(m_omega, weights, args.plot_output_dir)

    print(f"Writing weighted ROOT file to {output_root}")
    write_weighted_root(output_root, masses, weights)

    summary = {
        "skim_source": str(args.skim_source),
        "output_root": str(output_root),
        "mean": {"value": XIC_MEAN_PDG, "error": 0.0},
        "sigma": {"value": parameter_value(params["sigma"]), "error": parameter_error(fit_errors, params["sigma"])},
        "signal_yield": {"value": syield, "error": parameter_error(fit_errors, params["nsig"])},
        "background_yield": {"value": byield, "error": parameter_error(fit_errors, params["nbkg"])},
    }
    print(json.dumps(summary, indent=2))
    print("Done.")


if __name__ == "__main__":
    main()
