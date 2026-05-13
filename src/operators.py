import numpy as np


def random_orthogonal(n: int, rng: np.random.Generator) -> np.ndarray:
    """Create a deterministic random orthogonal matrix from a QR factorization."""
    q, r = np.linalg.qr(rng.standard_normal((n, n)))
    signs = np.sign(np.diag(r))
    signs[signs == 0] = 1
    return q * signs


def exponential_singular_values(n: int, alpha: float) -> np.ndarray:
    """Return singular values with condition number 10**alpha."""
    return 10.0 ** (-alpha * np.arange(n) / (n - 1))


def make_ill_conditioned_matrix(
    n: int,
    alpha: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Construct A = U diag(s) V^T with exponentially decaying singular values."""
    u = random_orthogonal(n, rng)
    v = random_orthogonal(n, rng)
    singular_values = exponential_singular_values(n, alpha)
    a = (u * singular_values) @ v.T
    return a, singular_values


def make_test_signal(n: int) -> np.ndarray:
    """Construct a 1D signal with low-frequency, high-frequency, and local features."""
    t = np.linspace(0.0, 1.0, n, endpoint=False)
    signal = np.sin(2 * np.pi * 3 * t)
    signal += 0.5 * np.sin(2 * np.pi * 17 * t)
    signal += 0.9 * ((t > 0.38) & (t < 0.48)).astype(float)
    signal -= 0.7 * ((t > 0.68) & (t < 0.73)).astype(float)
    return signal


def add_relative_noise(
    y: np.ndarray,
    rho: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Add noise with ||noise|| / ||y|| = rho."""
    noise = rng.standard_normal(y.shape)
    norm = np.linalg.norm(noise.ravel())
    if norm == 0:
        raise ValueError("Generated zero noise vector.")
    noise = noise / norm * (rho * np.linalg.norm(y.ravel()))
    return y + noise, noise
