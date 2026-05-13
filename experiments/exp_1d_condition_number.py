import json

import numpy as np
import pandas as pd

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.metrics import amplification_factor, relative_l2_error
from src.operators import add_relative_noise, make_ill_conditioned_matrix, make_test_signal
from src.plotting import plot_errorbar, plot_signal_recovery_with_zoom
from src.solvers import solve_pinv, solve_tikhonov_from_svd


def run_trial(n: int, alpha: float, rho: float, lam: float, seed: int) -> dict[str, float | np.ndarray]:
    rng = np.random.default_rng(seed)
    x = make_test_signal(n)
    a, _ = make_ill_conditioned_matrix(n, alpha, rng)
    b = a @ x
    b_tilde, noise = add_relative_noise(b, rho, rng)
    u, s, vt = np.linalg.svd(a, full_matrices=False)

    x_pinv = solve_pinv(a, b_tilde)
    x_tikh = solve_tikhonov_from_svd(u, s, vt, b_tilde, lam)

    return {
        "pinv_error": relative_l2_error(x_pinv, x),
        "tikhonov_error": relative_l2_error(x_tikh, x),
        "pinv_amplification": amplification_factor(x_pinv, x, noise, b),
        "tikhonov_amplification": amplification_factor(x_tikh, x, noise, b),
        "x": x,
        "x_pinv": x_pinv,
        "x_tikh": x_tikh,
    }


def main() -> None:
    ensure_output_dirs()
    n = 256
    rho = 1e-3
    lam = 1e-6
    alphas = [2, 4, 6, 8]
    repeats = 20

    rows = []
    example = None
    for alpha in alphas:
        for trial in range(repeats):
            result = run_trial(n, alpha, rho, lam, RNG_SEED + 1000 * int(alpha) + trial)
            rows.append(
                {
                    "alpha": alpha,
                    "condition_number": 10**alpha,
                    "trial": trial,
                    "rho": rho,
                    "lambda": lam,
                    "pinv_error": result["pinv_error"],
                    "tikhonov_error": result["tikhonov_error"],
                    "pinv_amplification": result["pinv_amplification"],
                    "tikhonov_amplification": result["tikhonov_amplification"],
                }
            )
            if alpha == 8 and trial == 0:
                example = result

    df = pd.DataFrame(rows)
    csv_path = RESULTS_DIR / "exp_1d_condition_number.csv"
    df.to_csv(csv_path, index=False)

    summary = (
        df.groupby("condition_number")[["pinv_error", "tikhonov_error"]]
        .agg(["mean", "std"])
        .reset_index()
    )
    kappas = summary["condition_number"].to_numpy()
    pinv_mean = summary[("pinv_error", "mean")].to_numpy()
    pinv_std = summary[("pinv_error", "std")].to_numpy()
    tikh_mean = summary[("tikhonov_error", "mean")].to_numpy()
    tikh_std = summary[("tikhonov_error", "std")].to_numpy()

    plot_errorbar(
        FIGURES_DIR / "signal" / "fig_1d_pinv_error_vs_kappa.png",
        "Pseudoinverse error grows with condition number",
        kappas,
        pinv_mean,
        pinv_std,
        "condition number kappa(A)",
        "relative L2 error",
        xscale="log",
    )
    plot_errorbar(
        FIGURES_DIR / "signal" / "fig_1d_tikhonov_error_vs_kappa.png",
        "Tikhonov error vs condition number",
        kappas,
        tikh_mean,
        tikh_std,
        "condition number kappa(A)",
        "relative L2 error",
        xscale="log",
    )

    if example is not None:
        plot_signal_recovery_with_zoom(
            FIGURES_DIR / "signal" / "fig_1d_signal_recovery_kappa_1e8.png",
            "1D recovery example, kappa(A)=1e8",
            example["x"],
            example["x_pinv"],
            example["x_tikh"],
        )

    metadata = {
        "seed": RNG_SEED,
        "n": n,
        "rho": rho,
        "lambda": lam,
        "alphas": alphas,
        "repeats": repeats,
        "csv": str(csv_path),
    }
    (RESULTS_DIR / "exp_1d_condition_number.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
