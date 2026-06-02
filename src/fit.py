import argparse
import json
import math
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import matplotlib
import numpy as np
import ROOT
from iminuit import Minuit

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
hep.style.use("LHCb2")
plt.rcParams.update(
    {
        "font.size": 11,
        "axes.labelsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 12,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 4.0,
        "ytick.major.size": 4.0,
        "xtick.minor.width": 0.35,
        "ytick.minor.width": 0.35,
        "xtick.minor.size": 2.5,
        "ytick.minor.size": 2.5,
        "lines.linewidth": 0.9,
    }
)


def add_lhcb_label(ax):
    ax.text(
        0.63,
        0.88,
        "LHCb Open Data",
        transform=ax.transAxes,
        fontsize=15,
        fontweight="bold",
        color="#1f3b73",
        ha="left",
        va="top",
        bbox=dict(boxstyle="round,pad=0.24", facecolor="white", edgecolor="#1f3b73", linewidth=1.0, alpha=0.96),
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Fit the weighted Omega_c spectrum.")
    parser.add_argument(
        "weighted_file",
        nargs="?",
        default=None,
        help="ROOT file produced by src/sweight.py.",
    )
    parser.add_argument(
        "--plot-output-dir",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "fit",
        help="Directory where the fit plot will be written if --plot-output-file is not set.",
    )
    parser.add_argument(
        "--plot-output-file",
        type=Path,
        default=None,
        help="Path where the fit plot PNG will be written.",
    )
    parser.add_argument(
        "--results-json",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "fit" / "omega_fit.json",
        help="Path where the fit summary JSON will be written.",
    )
    parser.add_argument(
        "--fit-log-file",
        type=Path,
        default=Path(__file__).resolve().parent.parent / "fit" / "omega_fit.minuit.log",
        help="Path where the fit log will be written.",
    )
    return parser.parse_args()


FIT_RANGE = (2970.0, 3500.0)
FIT_PLOT_BINS = 110
NUMBER_OF_SIGNALS = 7
SIGNAL_COLORS = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
]

SIGNAL_STARTS = [
    [2000, 3000, 20],
    [1500, 3050, 2],
    [1500, 3065, 7],
    [1500, 3090, 2],
    [2000, 3119, 4],
    [2000, 3185, 20],
    [200, 3327, 30],
]

SIGNAL_LIMITS = [
    [(1, np.inf), (2990, 3010), (5, 25)],
    [(1, np.inf), (3047, 3053), (1, 20)],
    [(1, np.inf), (3062, 3069), (1, 15)],
    [(1, np.inf), (3087, 3094), (1, 15)],
    [(1, np.inf), (3115, 3125), (1, 20)],
    [(1, np.inf), (3181, 3195), (3, 45)],
    [(1, np.inf), (3318, 3337), (5, 40)],
]

BACKGROUND_START = [10.0, 0.5, -0.002, 0.0, 2965.0]
BACKGROUND_LIMITS = [(1, np.inf), (0.0, 5.0), (-0.02, 0.0), (-1e-5, 1e-5), (2940.0, 2970.0)]


def load_weighted_data(weighted_file: Path) -> tuple[np.ndarray, np.ndarray]:
    df = ROOT.RDataFrame("DecayTree/DecayTree", str(weighted_file))
    data = df.AsNumpy(columns=["Omega_cst0_M", "sweight_signal"])
    m_omega = np.asarray(data["Omega_cst0_M"], dtype=np.float64)
    weights = np.asarray(data["sweight_signal"], dtype=np.float64)
    mask = (
        np.isfinite(m_omega)
        & np.isfinite(weights)
        & (m_omega >= FIT_RANGE[0])
        & (m_omega <= FIT_RANGE[1])
    )
    return m_omega[mask], weights[mask]


def threshold_background_shape(m: np.ndarray, a: float, b1: float, b2: float, m_thr: float) -> np.ndarray:
    lo, hi = FIT_RANGE
    x = np.asarray(m, dtype=np.float64)
    dm = x - m_thr
    dm_pos = np.clip(dm, 0.0, None)
    exponent = np.clip(b1 * dm + b2 * dm**2, -700, 700)
    return np.power(dm_pos, a) * np.exp(exponent) * np.heaviside(dm, 1.0)


def normalized_background_pdf(m: np.ndarray, a: float, b1: float, b2: float, m_thr: float) -> np.ndarray:
    lo, hi = FIT_RANGE
    raw = threshold_background_shape(m, a, b1, b2, m_thr)
    grid = np.linspace(lo, hi, 1500)
    norm = np.trapezoid(threshold_background_shape(grid, a, b1, b2, m_thr), grid)
    if norm <= 0 or not np.isfinite(norm):
        return np.zeros_like(raw)
    return raw / norm


def normalized_breit_wigner_pdf(m: np.ndarray, mean: float, gamma: float) -> np.ndarray:
    lo, hi = FIT_RANGE
    x = np.asarray(m, dtype=np.float64)
    half_gamma = 0.5 * gamma
    raw = (1.0 / np.pi) * half_gamma / ((x - mean) ** 2 + half_gamma**2)
    norm = (
        np.arctan2(2.0 * (hi - mean), gamma)
        - np.arctan2(2.0 * (lo - mean), gamma)
    ) / np.pi
    if norm <= 0 or not np.isfinite(norm):
        return np.zeros_like(raw)
    return raw / norm


def build_fit_model(m_omega: np.ndarray, weights: np.ndarray):
    total_weight = float(np.sum(weights))
    if not np.isfinite(total_weight) or total_weight <= 0:
        total_weight = float(len(m_omega))

    signal_weight_template = np.array([2000, 1500, 1500, 1500, 2000, 2000, 200], dtype=np.float64)
    signal_weight_template /= signal_weight_template.sum()
    signal_total = 0.55 * total_weight

    param_names = ["A0", "a", "b1", "b2", "m_thr"]
    for i in range(NUMBER_OF_SIGNALS):
        param_names.extend([f"A{i + 1}", f"mean{i + 1}", f"gamma{i + 1}"])

    start_values = list(BACKGROUND_START)
    for i, (_, mean_start, gamma_start) in enumerate(SIGNAL_STARTS):
        start_values.extend([float(signal_total * signal_weight_template[i]), mean_start, gamma_start])

    step_sizes = [
        2000.0,
        0.2,
        0.001,
        1.0e-6,
        1.0,
    ]
    for i, (_, _, gamma_start) in enumerate(SIGNAL_STARTS):
        step_sizes.extend([max(25.0, 0.15 * start_values[5 + 3 * i]), 0.5, max(1.0, 0.25 * gamma_start)])

    def model_density(m: np.ndarray, *params):
        a0, a, b1, b2, m_thr = params[:5]
        total = a0 * normalized_background_pdf(m, a, b1, b2, m_thr)
        for i in range(NUMBER_OF_SIGNALS):
            offset = 5 + 3 * i
            amp, mean, gamma = params[offset : offset + 3]
            total += amp * normalized_breit_wigner_pdf(m, mean, gamma)
        return total

    def nll(*params):
        density = model_density(m_omega, *params)
        if not np.all(np.isfinite(density)) or np.any(density <= 0):
            return np.inf
        total_yield = params[0] + sum(params[5 + 3 * i] for i in range(NUMBER_OF_SIGNALS))
        return total_yield - np.sum(weights * np.log(density))

    minuit = Minuit(nll, *start_values, name=param_names)
    minuit.errordef = 0.5
    minuit.tol = 1e-3
    minuit.print_level = 1
    minuit.strategy = 1
    for name, step in zip(param_names, step_sizes):
        minuit.errors[name] = step
    for name, lim in zip(param_names[:5], BACKGROUND_LIMITS):
        minuit.limits[name] = lim
    minuit.fixed["m_thr"] = True
    for i, peak_limits in enumerate(SIGNAL_LIMITS, start=1):
        for suffix, lim in zip(("A", "mean", "gamma"), peak_limits):
            minuit.limits[f"{suffix}{i}"] = lim

    return minuit, model_density


def save_fit_plot(m_omega, weights, minuit, model_density, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)
    hist, edges = np.histogram(m_omega, bins=FIT_PLOT_BINS, weights=weights, range=FIT_RANGE)
    hist_var, _ = np.histogram(m_omega, bins=FIT_PLOT_BINS, weights=np.square(weights), range=FIT_RANGE)
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_width = edges[1] - edges[0]
    errors = np.sqrt(np.maximum(hist_var, 1.0e-12))
    errors[errors == 0] = 1.0

    plot_x = np.linspace(FIT_RANGE[0], FIT_RANGE[1], 500)
    fitted_params = [value for value in minuit.values]
    plot_model = model_density(plot_x, *fitted_params) * bin_width
    plot_background = fitted_params[0] * normalized_background_pdf(plot_x, *fitted_params[1:5]) * bin_width

    fig, (ax_main, ax_pull) = plt.subplots(
        2,
        1,
        figsize=(9.4, 5.6),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.0},
    )

    ax_main.errorbar(
        centers,
        hist,
        errors,
        linestyle="",
        marker="x",
        ms=2.1,
        color="black",
        elinewidth=0.4,
        capsize=0,
        alpha=0.8,
        label="Data",
    )
    ax_main.plot(plot_x, plot_model, color="#1f3b73", linewidth=0.8, label="Total fit")
    ax_main.plot(plot_x, plot_background, linestyle="--", color="tab:green", linewidth=0.8, alpha=0.85, label="Background")
    signal_labels = [
        r"$\Omega_c^{**0}(3000)$",
        r"$\Omega_c^{**0}(3050)$",
        r"$\Omega_c^{**0}(3065)$",
        r"$\Omega_c^{**0}(3090)$",
        r"$\Omega_c^{**0}(3119)$",
        r"$\Omega_c^{**0}(3185)$",
        r"$\Omega_c^{**0}(3327)$",
    ]
    for i in range(NUMBER_OF_SIGNALS):
        offset = 5 + 3 * i
        amp, mean, gamma = fitted_params[offset : offset + 3]
        component = amp * normalized_breit_wigner_pdf(plot_x, mean, gamma) * bin_width
        ax_main.plot(
            plot_x,
            component,
            linestyle="-",
            color=SIGNAL_COLORS[i % len(SIGNAL_COLORS)],
            linewidth=1.3,
            alpha=0.95,
            label=signal_labels[i],
        )

    ax_main.set_ylabel(f"Weighted candidates / {bin_width:.2f} MeV")
    ax_main.set_xlim(FIT_RANGE)
    main_max = max(
        float(np.max(hist)),
        float(np.max(plot_model)),
        float(np.max(plot_background)),
    )
    ax_main.set_ylim(0, 1.15 * main_max)
    ax_main.grid(True, which="both", alpha=0.2, linewidth=0.5)
    ax_main.legend(fontsize=12, ncol=2, frameon=False, loc="lower right", handlelength=2.6, labelspacing=0.5)
    add_lhcb_label(ax_main)

    model_at_centers = model_density(centers, *fitted_params) * bin_width
    pull = (hist - model_at_centers) / errors
    ax_pull.bar(
        centers,
        pull,
        width=bin_width,
        align="center",
        color="black",
        edgecolor="black",
        linewidth=0.25,
    )
    ax_pull.axhline(0, color="gray", linestyle="-", linewidth=0.7, alpha=0.45)
    ax_pull.axhline(-3, color="red", linestyle="--", linewidth=0.8, alpha=0.45)
    ax_pull.axhline(3, color="red", linestyle="--", linewidth=0.8, alpha=0.45)
    ax_pull.set_ylabel("Pull [$\\sigma$]", labelpad=20)
    ax_pull.set_xlabel(r"$m(\Xi_c^+K^-)$ / MeV")
    ax_pull.set_ylim((-5, 5))
    ax_pull.set_xlim(FIT_RANGE)
    ax_pull.yaxis.tick_right()
    ax_pull.yaxis.set_label_position("left")
    ax_pull.tick_params(axis="y", labelright=True, labelleft=False)
    ax_pull.grid(True, which="both", alpha=0.2, linewidth=0.5)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def write_results_json(minuit: Minuit, output_path: Path, weighted_file: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary = {
        "weighted_file": str(weighted_file),
        "parameters": {},
    }
    for name in minuit.parameters:
        summary["parameters"][name] = {
            "value": float(minuit.values[name]),
            "error": float(minuit.errors[name]),
        }
    output_path.write_text(json.dumps(summary, indent=2))


def main():
    args = parse_args()
    weighted_file = (
        Path(args.weighted_file)
        if args.weighted_file is not None
        else Path(__file__).resolve().parent.parent / "weighted" / "omegac2xicK.root"
    )
    plot_output_file = args.plot_output_file or (args.plot_output_dir / "omega_fit.png")

    print(f"Reading weighted file {weighted_file}")
    m_omega, weights = load_weighted_data(weighted_file)
    print(f"Fitting weighted Omega_c spectrum with {len(m_omega)} entries")

    log_path = args.fit_log_file
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log_file, redirect_stdout(log_file), redirect_stderr(log_file):
        minuit, model_density = build_fit_model(m_omega, weights)
        minuit.simplex(ncall=10000)
        minuit.strategy = 2
        minuit.migrad(ncall=100000)
        if not minuit.fmin.is_valid:
            minuit.migrad(ncall=100000)
        minuit.hesse()
        print(minuit)
    print(f"Wrote fit log to {log_path}")

    print("Writing fit plot")
    save_fit_plot(m_omega, weights, minuit, model_density, plot_output_file)
    write_results_json(minuit, args.results_json, weighted_file)
    print(f"Wrote fit summary to {args.results_json}")
    print(f"Wrote fit plot to {plot_output_file}")
    print("Done.")


if __name__ == "__main__":
    main()
