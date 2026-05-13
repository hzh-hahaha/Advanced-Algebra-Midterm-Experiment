import numpy as np


def relative_l2_error(x_hat: np.ndarray, x_true: np.ndarray) -> float:
    return float(np.linalg.norm((x_hat - x_true).ravel()) / np.linalg.norm(x_true.ravel()))


def relative_residual(a: np.ndarray, x_hat: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm((a @ x_hat - b).ravel()) / np.linalg.norm(b.ravel()))


def frobenius_relative_error(x_hat: np.ndarray, x_true: np.ndarray) -> float:
    return float(np.linalg.norm(x_hat - x_true, ord="fro") / np.linalg.norm(x_true, ord="fro"))


def amplification_factor(x_hat: np.ndarray, x_true: np.ndarray, noise: np.ndarray, b: np.ndarray) -> float:
    output_error = np.linalg.norm((x_hat - x_true).ravel()) / np.linalg.norm(x_true.ravel())
    input_error = np.linalg.norm(noise.ravel()) / np.linalg.norm(b.ravel())
    return float(output_error / input_error)
