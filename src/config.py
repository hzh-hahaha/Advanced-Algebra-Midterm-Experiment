from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"
DATA_DIR = PROJECT_ROOT / "data"

RNG_SEED = 20260513


def ensure_output_dirs() -> None:
    for path in [
        FIGURES_DIR / "theory",
        FIGURES_DIR / "signal",
        FIGURES_DIR / "image",
        RESULTS_DIR,
        DATA_DIR / "raw",
        DATA_DIR / "processed",
    ]:
        path.mkdir(parents=True, exist_ok=True)
