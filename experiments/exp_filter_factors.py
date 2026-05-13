import json

import matplotlib.pyplot as plt
import numpy as np

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.operators import exponential_singular_values
from src.solvers import tikhonov_filter_factors


def main() -> None:
    ensure_output_dirs()
    n = 256
    alpha = 8
    lambdas = np.array([1e-12, 1e-8, 1e-4, 1e-2])
    s = exponential_singular_values(n, alpha)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    axes[0].plot(np.arange(1, n + 1), s, linewidth=1.8)
    axes[0].set_yscale("log")
    axes[0].set_title("Singular value spectrum")
    axes[0].set_xlabel("index i")
    axes[0].set_ylabel("sigma_i")
    axes[0].grid(alpha=0.25, which="both")

    axes[1].plot(np.arange(1, n + 1), 1 / s, label="pseudoinverse 1/sigma_i", linewidth=1.8)
    for lam in lambdas:
        axes[1].plot(
            np.arange(1, n + 1),
            tikhonov_filter_factors(s, lam),
            label=f"Tikhonov lambda={lam:g}",
            linewidth=1.5,
        )
    axes[1].set_yscale("log")
    axes[1].set_title("Spectral amplification factors")
    axes[1].set_xlabel("index i")
    axes[1].set_ylabel("amplification")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25, which="both")

    fig.tight_layout()
    out = FIGURES_DIR / "theory" / "fig_filter_factors.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)

    metadata = {
        "seed": RNG_SEED,
        "n": n,
        "alpha": alpha,
        "condition_number": float(10**alpha),
        "lambdas": lambdas.tolist(),
        "figure": str(out),
    }
    (RESULTS_DIR / "exp_filter_factors.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
