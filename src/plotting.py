from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def savefig(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path, dpi=220)
    plt.close()


def plot_signal_comparison(path: Path, title: str, series: dict[str, np.ndarray]) -> None:
    plt.figure(figsize=(9, 4.8))
    for label, values in series.items():
        plt.plot(values, label=label, linewidth=1.6)
    plt.title(title)
    plt.xlabel("index")
    plt.ylabel("value")
    plt.legend()
    plt.grid(alpha=0.25)
    savefig(path)


def plot_signal_recovery_with_zoom(
    path: Path,
    title: str,
    x_true: np.ndarray,
    x_pinv: np.ndarray,
    x_tikh: np.ndarray,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(x_true, label="true x", linewidth=1.8)
    axes[0].plot(x_pinv, label="pseudoinverse", linewidth=1.3)
    axes[0].plot(x_tikh, label="Tikhonov", linewidth=1.8)
    axes[0].set_title("Full scale")
    axes[0].set_xlabel("index")
    axes[0].set_ylabel("value")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.25)

    axes[1].plot(x_true, label="true x", linewidth=1.8)
    axes[1].plot(x_tikh, label="Tikhonov", linewidth=1.8)
    ymax = max(float(np.max(np.abs(x_true))), float(np.max(np.abs(x_tikh))))
    axes[1].set_ylim(-1.25 * ymax, 1.25 * ymax)
    axes[1].set_title("Zoomed view without pseudoinverse blow-up")
    axes[1].set_xlabel("index")
    axes[1].set_ylabel("value")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.25)

    fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_errorbar(path: Path, title: str, x, y, yerr, xlabel: str, ylabel: str, xscale: str = "linear") -> None:
    plt.figure(figsize=(7, 4.8))
    plt.errorbar(x, y, yerr=yerr, marker="o", capsize=4, linewidth=1.7)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale(xscale)
    plt.yscale("log")
    plt.grid(alpha=0.25, which="both")
    savefig(path)


def plot_lambda_curve(path: Path, lambdas: np.ndarray, errors: np.ndarray, residuals: np.ndarray, norms: np.ndarray) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    axes[0].plot(lambdas, errors, marker=".", linewidth=1.4)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_title("Relative error vs lambda")
    axes[0].set_xlabel("lambda")
    axes[0].set_ylabel("relative error")
    axes[0].grid(alpha=0.25, which="both")

    axes[1].plot(residuals, norms, marker=".", linewidth=1.4)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_title("L-curve")
    axes[1].set_xlabel("||A x_lambda - b_tilde||_2")
    axes[1].set_ylabel("||x_lambda||_2")
    axes[1].grid(alpha=0.25, which="both")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_image_grid(path: Path, title: str, images: list[tuple[str, np.ndarray]]) -> None:
    cols = len(images)
    fig, axes = plt.subplots(1, cols, figsize=(3.2 * cols, 3.5))
    if cols == 1:
        axes = [axes]
    for ax, (label, image) in zip(axes, images):
        ax.imshow(image, cmap="gray")
        ax.set_title(label)
        ax.axis("off")
    fig.suptitle(title)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)
