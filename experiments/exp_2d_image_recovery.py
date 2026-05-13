import json

import numpy as np
import pandas as pd
from skimage import data, img_as_float
from skimage.metrics import peak_signal_noise_ratio, structural_similarity
from skimage.transform import resize

from src.config import FIGURES_DIR, RESULTS_DIR, RNG_SEED, ensure_output_dirs
from src.metrics import frobenius_relative_error
from src.operators import add_relative_noise, make_ill_conditioned_matrix
from src.plotting import plot_image_grid


def normalize_for_display(x: np.ndarray) -> np.ndarray:
    lo = np.percentile(x, 1)
    hi = np.percentile(x, 99)
    if hi <= lo:
        return np.clip(x, 0, 1)
    return np.clip((x - lo) / (hi - lo), 0, 1)


def tikhonov_image_recover(
    b_mat: np.ndarray,
    y_tilde: np.ndarray,
    lam: float,
) -> np.ndarray:
    a = np.kron(b_mat, b_mat)
    y = y_tilde.reshape(-1)
    lhs = a.T @ a + lam * np.eye(a.shape[1])
    rhs = a.T @ y
    return np.linalg.solve(lhs, rhs).reshape(y_tilde.shape)


def main() -> None:
    ensure_output_dirs()
    m = 64
    alpha_b = 3
    rho = 1e-3
    lambdas = [1e-10, 1e-7, 1e-5, 1e-3]
    rng = np.random.default_rng(RNG_SEED + 222)

    image = img_as_float(data.camera())
    x = resize(image, (m, m), anti_aliasing=True)
    b_mat, _ = make_ill_conditioned_matrix(m, alpha_b, rng)
    y = b_mat @ x @ b_mat.T
    y_tilde, noise = add_relative_noise(y, rho, rng)

    x_pinv = np.linalg.pinv(b_mat, rcond=1e-15) @ y_tilde @ np.linalg.pinv(b_mat, rcond=1e-15).T

    rows = [
        {
            "method": "pseudoinverse",
            "lambda": np.nan,
            "frobenius_relative_error": frobenius_relative_error(x_pinv, x),
            "psnr": peak_signal_noise_ratio(x, np.clip(x_pinv, 0, 1), data_range=1.0),
            "ssim": structural_similarity(x, np.clip(x_pinv, 0, 1), data_range=1.0),
        }
    ]

    recovered = []
    for lam in lambdas:
        x_lam = tikhonov_image_recover(b_mat, y_tilde, lam)
        recovered.append((lam, x_lam))
        rows.append(
            {
                "method": "tikhonov",
                "lambda": lam,
                "frobenius_relative_error": frobenius_relative_error(x_lam, x),
                "psnr": peak_signal_noise_ratio(x, np.clip(x_lam, 0, 1), data_range=1.0),
                "ssim": structural_similarity(x, np.clip(x_lam, 0, 1), data_range=1.0),
            }
        )

    df = pd.DataFrame(rows)
    csv_path = RESULTS_DIR / "exp_2d_image_recovery.csv"
    df.to_csv(csv_path, index=False)

    best = df[df["method"] == "tikhonov"].sort_values("frobenius_relative_error").iloc[0]
    best_image = next(img for lam, img in recovered if lam == best["lambda"])

    plot_image_grid(
        FIGURES_DIR / "image" / "fig_2d_image_recovery_main.png",
        "2D ill-conditioned image recovery",
        [
            ("true X", x),
            ("degraded BXB^T", normalize_for_display(y)),
            ("noisy observation", normalize_for_display(y_tilde)),
            ("pseudoinverse", normalize_for_display(x_pinv)),
            (f"Tikhonov lambda={best['lambda']:.0e}", normalize_for_display(best_image)),
        ],
    )

    plot_image_grid(
        FIGURES_DIR / "image" / "fig_2d_lambda_grid.png",
        "Tikhonov recovery under different lambda values",
        [(f"lambda={lam:.0e}", normalize_for_display(img)) for lam, img in recovered],
    )

    metadata = {
        "seed": RNG_SEED + 222,
        "image_size": m,
        "alpha_B": alpha_b,
        "condition_number_B": 10**alpha_b,
        "condition_number_kron_B_B": 10 ** (2 * alpha_b),
        "rho": rho,
        "noise_norm": float(np.linalg.norm(noise, ord="fro")),
        "lambdas": lambdas,
        "best_lambda": float(best["lambda"]),
        "best_frobenius_relative_error": float(best["frobenius_relative_error"]),
        "csv": str(csv_path),
    }
    (RESULTS_DIR / "exp_2d_image_recovery.json").write_text(
        json.dumps(metadata, indent=2),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
