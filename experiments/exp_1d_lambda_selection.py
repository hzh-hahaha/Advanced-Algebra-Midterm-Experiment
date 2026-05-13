import json

import numpy as np
import pandas as pd

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.metrics import relative_l2_error
from src.operators import add_relative_noise, make_ill_conditioned_matrix, make_test_signal
from src.plotting import plot_lambda_curve
from src.solvers import solve_tikhonov_from_svd


def main() -> None:
    ensure_output_dirs()
    n = 256
    alpha = 8
    rho = 1e-3
    lambdas = np.logspace(-16, 0, 120)

    rng = np.random.default_rng(RNG_SEED + 77)
    x = make_test_signal(n)
    a, _ = make_ill_conditioned_matrix(n, alpha, rng)
    b = a @ x
    b_tilde, noise = add_relative_noise(b, rho, rng)
    u, s, vt = np.linalg.svd(a, full_matrices=False)

    rows = []
    for lam in lambdas:
        x_hat = solve_tikhonov_from_svd(u, s, vt, b_tilde, lam)
        rows.append(
            {
                "lambda": lam,
                "relative_error": relative_l2_error(x_hat, x),
                "residual_norm": float(np.linalg.norm(a @ x_hat - b_tilde)),
                "solution_norm": float(np.linalg.norm(x_hat)),
                "noise_norm": float(np.linalg.norm(noise)),
            }
        )

    df = pd.DataFrame(rows)
    best = df.loc[df["relative_error"].idxmin()]
    csv_path = RESULTS_DIR / "exp_1d_lambda_selection.csv"
    df.to_csv(csv_path, index=False)

    plot_lambda_curve(
        FIGURES_DIR / "signal" / "fig_1d_lambda_selection_and_lcurve.png",
        df["lambda"].to_numpy(),
        df["relative_error"].to_numpy(),
        df["residual_norm"].to_numpy(),
        df["solution_norm"].to_numpy(),
    )

    metadata = {
        "seed": RNG_SEED + 77,
        "n": n,
        "alpha": alpha,
        "condition_number": 10**alpha,
        "rho": rho,
        "lambda_min": float(lambdas.min()),
        "lambda_max": float(lambdas.max()),
        "lambda_count": len(lambdas),
        "oracle_lambda": float(best["lambda"]),
        "oracle_relative_error": float(best["relative_error"]),
        "noise_norm": float(np.linalg.norm(noise)),
        "csv": str(csv_path),
    }
    (RESULTS_DIR / "exp_1d_lambda_selection.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
