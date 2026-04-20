"""
Generate assets for the Wolfram community post and README:
- Download real Horsehead Nebula FITS cutouts from NASA SkyView.
- Render percentile-clipped PNG previews.
- Generate a simulated Horsehead image (always available, network-free).
- Produce a 4-panel analysis figure.
"""

import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
from io import BytesIO

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

HORSEHEAD_RA = 85.2479
HORSEHEAD_DEC = -2.4583


def try_download(survey, fov=0.5, pixels=600):
    try:
        import requests
        from astropy.io import fits as afits

        url = "https://skyview.gsfc.nasa.gov/cgi-bin/images"
        params = {
            "Position": f"{HORSEHEAD_RA},{HORSEHEAD_DEC}",
            "Survey":   survey,
            "Pixels":   f"{pixels},{pixels}",
            "Size":     fov,
            "Return":   "FITS",
            "Scaling":  "Linear",
        }
        r = requests.get(url, params=params, timeout=45)
        if r.status_code == 200 and len(r.content) > 2000:
            with afits.open(BytesIO(r.content)) as hdul:
                data = hdul[0].data.astype(float)
            return r.content, data
    except Exception as e:  # noqa: BLE001
        print(f"  ! {survey} download failed: {e}")
    return None, None


def save_preview(data, path, title, cmap="gray"):
    vmin, vmax = np.nanpercentile(data, [1, 99.5])
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.imshow(data, cmap=cmap, origin="lower", vmin=vmin, vmax=vmax,
              interpolation="nearest")
    ax.set_title(title, fontsize=12)
    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")
    plt.tight_layout()
    fig.savefig(path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {path}")


def simulate():
    size = 600
    x = np.linspace(-1, 1, size)
    y = np.linspace(-1, 1, size)
    X, Y = np.meshgrid(x, y)
    rng = np.random.default_rng(42)

    bg = 100 + 150 * np.exp(-(X ** 2 + Y ** 2) / 0.5)
    for _ in range(5):
        a = rng.uniform(0, np.pi)
        f = rng.uniform(5, 15)
        A = rng.uniform(10, 30)
        bg += A * np.sin(f * (X * np.cos(a) + Y * np.sin(a)))

    mask = np.exp(-(((X - 0.10) / 0.30) ** 2 +
                    ((Y + 0.20) / 0.40) ** 2) * 8)
    mask += 0.7 * np.exp(-(((X - 0.15) / 0.15) ** 2 +
                           ((Y - 0.15) / 0.35) ** 2) * 6)
    tau = 3.0 * mask
    img = bg * np.exp(-tau)
    img = rng.poisson(np.clip(img, 0, None)).astype(float)
    img += rng.normal(0, 0.02 * img.mean(), img.shape)
    img = np.maximum(img, 0)
    return img, bg, tau


def analysis_panel(img, out_path):
    threshold = np.percentile(img, 30)
    dark = img < threshold

    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    vmin, vmax = np.nanpercentile(img, [1, 99.5])
    axes[0, 0].imshow(img, cmap="gray", origin="lower",
                      vmin=vmin, vmax=vmax)
    axes[0, 0].set_title("Raw image")
    axes[0, 0].axis("off")

    axes[0, 1].imshow(np.log10(img + 1), cmap="inferno", origin="lower")
    axes[0, 1].set_title("log10 scaling")
    axes[0, 1].axis("off")

    masked = np.ma.masked_where(~dark, img)
    axes[1, 0].imshow(img, cmap="gray", origin="lower",
                      vmin=vmin, vmax=vmax, alpha=0.6)
    axes[1, 0].imshow(masked, cmap="Reds", origin="lower", alpha=0.9)
    axes[1, 0].set_title("Dark-nebula segmentation (bottom 30%)")
    axes[1, 0].axis("off")

    axes[1, 1].hist(img.flatten(), bins=80, color="steelblue",
                    edgecolor="black")
    axes[1, 1].axvline(threshold, color="red", ls="--",
                       label=f"30th pct = {threshold:.1f}")
    axes[1, 1].set_xlabel("Pixel intensity (ADU)")
    axes[1, 1].set_ylabel("Count")
    axes[1, 1].set_title("Intensity distribution")
    axes[1, 1].legend()

    fig.suptitle("Horsehead Nebula - rudimentary analysis", fontsize=14)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    os.chdir(OUT_DIR)
    print("Generating community-post assets in", OUT_DIR)

    any_real = False
    real_data = None
    for survey in ("DSS2 Red", "2MASS-J"):
        print(f"- trying SkyView survey {survey!r}")
        content, data = try_download(survey, fov=0.5, pixels=600)
        if content is not None:
            safe = survey.replace(" ", "_").replace("-", "_")
            fits_path = os.path.join(OUT_DIR, f"horsehead_{safe}.fits")
            png_path = os.path.join(OUT_DIR, f"horsehead_{safe}.png")
            with open(fits_path, "wb") as f:
                f.write(content)
            save_preview(data, png_path, f"Horsehead Nebula - {survey}")
            if real_data is None:
                real_data = data
            any_real = True

    print("- generating simulated Horsehead image")
    sim_img, bg, tau = simulate()
    np.save(os.path.join(OUT_DIR, "horsehead_simulated.npy"), sim_img)
    save_preview(sim_img, os.path.join(OUT_DIR, "horsehead_simulated.png"),
                 "Simulated Horsehead Nebula", cmap="gray")

    analysis_src = real_data if real_data is not None else sim_img
    analysis_panel(analysis_src,
                   os.path.join(OUT_DIR, "horsehead_analysis_panel.png"))

    # Dump stats as JSON for the post to reference
    threshold = float(np.percentile(analysis_src, 30))
    dark = analysis_src < threshold
    bg_mean = float(np.mean(analysis_src[~dark]))
    neb_mean = float(np.mean(analysis_src[dark]))
    tau_est = float(-np.log((neb_mean + 1e-6) / (bg_mean + 1e-6)))
    stats = {
        "source": "real" if any_real else "simulated",
        "shape": list(analysis_src.shape),
        "min": float(np.min(analysis_src)),
        "max": float(np.max(analysis_src)),
        "mean": float(np.mean(analysis_src)),
        "median": float(np.median(analysis_src)),
        "threshold_p30": threshold,
        "background_mean": bg_mean,
        "dark_mean": neb_mean,
        "estimated_optical_depth": tau_est,
        "estimated_extinction_mag": tau_est * 1.086,
    }
    with open(os.path.join(OUT_DIR, "analysis_stats.json"), "w") as f:
        json.dump(stats, f, indent=2)
    print("  wrote analysis_stats.json")

    print("done; real data:", any_real)


if __name__ == "__main__":
    main()
