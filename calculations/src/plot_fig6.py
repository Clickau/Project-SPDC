"""
plot_fig6.py -- recreate Fig. 6(e)-(h) of Hsu & Lai (arXiv:1302.6226).

Four stacked log-scale color maps of the angle-resolved SPDC flux
Ns(lam_s, theta_s) in BBO, one per phase-matching angle theta_m, with the
perfect-phase-matching tuning curve overlaid as a white dashed line
(as in the paper).  theta_s is the EXTERNAL (air) angle, matching the
paper's y-axis.

All physics (corrected Eq. 9 flux, Sellmeier, tuning curve) lives in
spdc_physics.py.  Writes spdc_fig6.png.

Note: with our Sellmeier coefficients the collinear-degenerate point
(810 nm at theta_s = 0) lands at theta_m ~ 28.82 deg vs the paper's 28.6 deg
(the paper used a NIST index program) -- a ~0.2 deg label offset; the physics
is unchanged.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

import spdc_physics as phys

OUTDIR = Path(__file__).resolve().parent.parent / "output"


def make_figure(lam_min=430e-9, lam_max=1000e-9, dlam=1e-9,
                th_max_deg=6.0, dth_mrad=0.5,
                thetams_deg=(28.6, 28.8, 29.1, 29.4),
                fname=OUTDIR / "spdc_fig6.png"):
    lam_grid = np.arange(lam_min, lam_max + 1e-12, dlam)             # m
    dth = dth_mrad * 1e-3                                            # rad
    th_grid = np.arange(-th_max_deg, th_max_deg + 1e-9,
                        np.degrees(dth)) * np.pi / 180               # rad

    lam_nm = lam_grid * 1e9
    th_deg = np.degrees(th_grid)

    fig, axes = plt.subplots(len(thetams_deg), 1, figsize=(7, 11), sharex=True)

    for ax, tm_deg in zip(axes, thetams_deg):
        theta_m = np.radians(tm_deg)

        # Build the map column by column (one wavelength at a time) to keep
        # the memory of the xi_i quadrature axis bounded.
        img = np.zeros((th_grid.size, lam_grid.size))
        for j, lam_s in enumerate(lam_grid):
            img[:, j] = phys.Ns(lam_s, th_grid, theta_m, dlam, dth)

        img = np.maximum(img, 1.0)   # 1 count/s floor so the log scale behaves
        pcm = ax.pcolormesh(lam_nm, th_deg, img,
                            norm=LogNorm(vmin=1, vmax=1e6),
                            cmap="jet", shading="auto")

        # perfect phase-matching tuning curve (white dashed, like the paper)
        tc = phys.tuning_curve(lam_grid, theta_m)
        ax.plot(lam_nm,  tc, "w--", lw=0.8)
        ax.plot(lam_nm, -tc, "w--", lw=0.8)

        ax.set_ylabel(r"$\theta_s$ (deg)")
        ax.set_ylim(-th_max_deg, th_max_deg)
        ax.text(0.02, 0.86, rf"$\theta_m = {tm_deg:.1f}^\circ$",
                transform=ax.transAxes, color="white", fontsize=11)
        fig.colorbar(pcm, ax=ax, pad=0.01, label="counts/s")

    axes[-1].set_xlabel(r"$\lambda_s$ (nm)")
    fig.suptitle("Recreated Fig. 6(e)-(h): SPDC flux $N_s(\\lambda_s,\\theta_s)$ in BBO",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(fname, dpi=140)
    print(f"saved {fname}")


if __name__ == "__main__":
    phys.sanity_check()
    make_figure()
