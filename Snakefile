configfile: "config.yaml"

import sys
from pathlib import Path
from snakemake.io import directory

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
from skim import collect_files

def parse_bool(value, default=True):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


RUN_SKIM = parse_bool(config.get("run_skim", config["workflow"].get("run_skim", True)))
RAW_SOURCE = config["input"]["raw_source"]
MAX_FILES = int(config.get("max_files", config["input"].get("max_files", 0)))
MAX_FILES_LABEL = "ALL" if MAX_FILES <= 0 else str(MAX_FILES)
SKIMMED_DIR = str(Path(config["paths"]["skimmed_dir"]))
EOS_SKIMMED_DIR = config["paths"]["eos_skimmed_dir"]
PLOTS_DIR = str(Path(config["paths"]["plots_dir"]))
WEIGHTED_DIR = str(Path(config["paths"]["weighted_dir"]))
FIT_DIR = str(Path(config["paths"]["fit_dir"]))
PLOTS_MAX_DIR = str(Path(PLOTS_DIR) / f"maxfiles_{MAX_FILES_LABEL}")
WEIGHTED_MAX_DIR = str(Path(WEIGHTED_DIR) / f"maxfiles_{MAX_FILES_LABEL}")
FIT_MAX_DIR = str(Path(FIT_DIR) / f"maxfiles_{MAX_FILES_LABEL}")


def skimmed_file() -> str:
    return SKIMMED_DIR


def weighted_file() -> str:
    return str(Path(WEIGHTED_MAX_DIR) / "omegac2xicK.root")


def plot_dir(name: str) -> str:
    return str(Path(PLOTS_MAX_DIR) / name)


def fit_result_file() -> str:
    return str(Path(FIT_MAX_DIR) / "omega_fit.json")


def fit_plot_file() -> str:
    return str(Path(FIT_MAX_DIR) / "omega_fit.png")


def fit_log_file() -> str:
    return str(Path(FIT_MAX_DIR) / "omega_fit.minuit.log")


def skim_dir(_wildcards):
    if RUN_SKIM:
        return checkpoints.skim.get().output[0]
    return EOS_SKIMMED_DIR


ALL_TARGETS = [
    plot_dir("skim/xic_mass.png"),
    plot_dir("skim/omega_c_mass.png"),
    weighted_file(),
    fit_result_file(),
    fit_plot_file(),
]
if RUN_SKIM:
    ALL_TARGETS.append(directory(skimmed_file()))


rule all:
    input:
        ALL_TARGETS,


checkpoint skim:
    output:
        directory(skimmed_file()),
    params:
        source=RAW_SOURCE,
        max_files=MAX_FILES,
        skimmed_dir=SKIMMED_DIR,
    shell:
        "python3 src/skim.py {params.source:q} --max-files {params.max_files} --output-dir {params.skimmed_dir:q}"


if RUN_SKIM:

    rule plot:
        input:
            skim_dir=skim_dir,
        output:
            plot_dir("skim/xic_mass.png"),
            plot_dir("skim/omega_c_mass.png"),
        params:
            plots_dir=PLOTS_MAX_DIR,
            max_files=MAX_FILES,
        shell:
            "python3 src/skim_plots.py {input.skim_dir:q} --max-files {params.max_files} --plot-output-dir {params.plots_dir:q}/skim"


    rule sweight:
        input:
            skim_dir=skim_dir,
        output:
            weighted_file(),
        params:
            plots_dir=PLOTS_MAX_DIR,
            max_files=MAX_FILES,
        shell:
            "python3 src/sweight.py {input.skim_dir:q} --max-files {params.max_files} --output-root {output[0]:q} --plot-output-dir {params.plots_dir:q}/sweight"

else:

    rule plot:
        output:
            plot_dir("skim/xic_mass.png"),
            plot_dir("skim/omega_c_mass.png"),
        params:
            skim_source=EOS_SKIMMED_DIR,
            plots_dir=PLOTS_MAX_DIR,
            max_files=MAX_FILES,
        shell:
            "python3 src/skim_plots.py {params.skim_source:q} --max-files {params.max_files} --plot-output-dir {params.plots_dir:q}/skim"


    rule sweight:
        output:
            weighted_file(),
        params:
            skim_source=EOS_SKIMMED_DIR,
            plots_dir=PLOTS_MAX_DIR,
            max_files=MAX_FILES,
        shell:
            "python3 src/sweight.py {params.skim_source:q} --max-files {params.max_files} --output-root {output[0]:q} --plot-output-dir {params.plots_dir:q}/sweight"


rule fit:
    input:
        weighted_file(),
    output:
        fit_result_file(),
        fit_plot_file(),
    params:
        fit_log=fit_log_file(),
    shell:
        "python3 src/fit.py {input[0]:q} --results-json {output[0]:q} --plot-output-file {output[1]:q} --fit-log-file {params.fit_log:q}"
