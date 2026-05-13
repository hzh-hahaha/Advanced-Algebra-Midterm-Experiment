import numpy as np


def solve_pinv(a: np.ndarray, b: np.ndarray, rcond: float = 1e-15) -> np.ndarray:
    return np.linalg.pinv(a, rcond=rcond) @ b


def solve_tikhonov(a: np.ndarray, b: np.ndarray, lam: float) -> np.ndarray:
    n = a.shape[1]
    lhs = a.T @ a + lam * np.eye(n)
    rhs = a.T @ b
    return np.linalg.solve(lhs, rhs)


def solve_tikhonov_from_svd(
    u: np.ndarray,
    s: np.ndarray,
    vt: np.ndarray,
    b: np.ndarray,
    lam: float,
) -> np.ndarray:
    coeff = (s / (s * s + lam)) * (u.T @ b)
    return vt.T @ coeff


def solve_truncated_svd(
    u: np.ndarray,
    s: np.ndarray,
    vt: np.ndarray,
    b: np.ndarray,
    rank: int,
) -> np.ndarray:
    coeff = (u[:, :rank].T @ b) / s[:rank]
    return vt[:rank, :].T @ coeff


def tikhonov_filter_factors(s: np.ndarray, lam: float) -> np.ndarray:
    return s / (s * s + lam)
