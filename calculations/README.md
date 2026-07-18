# SPDC-in-BBO flux calculations

Numerical model of the angle-resolved parametric-fluorescence (SPDC) photon
flux from a BBO crystal, built on **Hsu & Lai, "Angle-Resolved Spectroscopy of
Parametric Fluorescence" ([arXiv:1302.6226](https://arxiv.org/abs/1302.6226))**.
Type-I phase matching (e → o + o), pump λp = 405 nm, degenerate signal/idler
at 810 nm.

The code does two things:

1. **Reproduces the paper** — the angle-resolved flux maps of Fig. 6 and the
   paper's Table I efficiency numbers, from the paper's Eq. (9).
2. **Predicts absolute count rates for a real tabletop setup** — singles on a
   finite-aperture detector and signal–idler coincidences between two
   detectors at ±3°, including measured/datasheet component efficiencies
   (mirrors, filter, fiber coupling, APD).

## Quick start

```
python src/spdc_physics.py        # sanity check against the paper (fast)
python src/plot_fig6.py           # recreate Fig. 6(e)-(h)      -> output/spdc_fig6*.png
python src/detector_counts.py     # singles rate, my setup      -> output/detector_scan_mysetup.png
python src/coincidence_counts.py  # coincidences, my setup      -> output/coincidence_scan.png
```

Requires Python 3 with `numpy` and `matplotlib`. Scripts can be run from the
repo root or from inside `src/` — output paths are resolved relative to the
script file. Every run also tees its printed summary (rates, peak angles,
sanity numbers) into a `.txt` next to the figure in `output/`.

## Layout

```
src/                    all source (physics library + runnable scripts)
src/archive/            superseded/testing scripts, LaTeX note, original Mathematica code
output/                 generated figures + teed text logs
output/archive/         outputs of the archived scripts
paper.pdf               the reference paper (Hsu & Lai)
mmf_focusing_findings.md  memo: pump-focusing experiment vs. model (future work)
CLAUDE.md               working notes / project instructions
```

---

## The main scripts

### `src/spdc_physics.py` — the physics library

All physics lives here; the runnable scripts contain none of their own.
Provides the SI constants, the paper's experiment parameters, the BBO
Sellmeier indices, and:

- **`Ns(lam_s, theta_s, theta_m, dlam, dtheta, ...)`** — the master equation:
  signal photon flux (counts/s) into one (dλ × dθ) bin, integrated over the
  full 2π azimuth of the emission cone. This one function serves the Fig. 6
  maps, the detector scans, and (via `idler_theta_window`) the coincidence
  calculation.
- **`tuning_curve(lam_array, theta_m)`** — external emission angle of exact
  phase matching vs signal wavelength (the white dashed curves in Fig. 6).
- **`crystal_rotation(theta_m, theta_m_cut)`** — external crystal-rotation
  angle α needed to reach an internal phase-matching angle θm, accounting for
  pump refraction at the entrance face.
- **`sanity_check()`** — run `python src/spdc_physics.py` to verify against
  the paper: η = 2Ns/Np ≈ 2.7e−10 and total degenerate flux 2Ns ≈ 4.4e7 /s at
  80 mW (paper: ~2.6–2.7e−10 and ~4.2e7 /s).

#### The master equation (paper Eq. 9)

`Ns` evaluates the paper's Appendix intermediate equation, changes variables
ω→λ and κ→θ, and integrates the signal azimuth over 2π. In the paper's own
variables it equals exactly **2π × the printed Eq. (9)**:

```
Ns(λs,θs) = 2π · [ 2π·√(2π) · ħ·deff²·ωp·L²·Np / (ε0·ns·ni·np·λs⁵·λi)
            · sin(2θs) · dλs · dθs
            · ∫dξi  exp(−½(ξs−ξi)²) · sinc²(½·L·Δkz) ]
```

with `sinc(x) = sin(x)/x`, the longitudinal phase mismatch
`Δkz = √(ks² − κs²) + √(ki² − (ξi/W)²) − kp` (`k = 2πn/λ`, in-medium),
dimensionless transverse momenta ξ = κ·W, and the idler fixed by energy
conservation, `λi = λs·λp/(λs − λp)`. The Gaussian `exp(−½(ξs−ξi)²)` is the
transverse phase-matching function of a Gaussian pump of 1/e² radius W.

Implementation details:

- **Angle convention:** θs is the **external (lab) angle** — κs = k0·sin θs
  with k0 = 2π/λs, the paper's Fig. 2 / Eq. (9) / Fig. 6 convention
  (`external=True`, the default). `external=False` parametrises the same κs
  by the internal angle; both give identical total flux (checked in
  `sanity_check()`), only the per-angle y-scale differs by ~ns.
- **ξi quadrature:** the Gaussian confines the integrand to ξi ≈ ξs, so the
  quadrature runs on a ±8σ window tracking ξs (`Nxi=200` trapezoid points by
  default). Evanescent idlers (κi > ki) are zeroed.
- **`lam_i_max`** — optional idler-wavelength cutoff (used at 3500 nm =
  BBO transparency edge / limit of the index data). Without it, short-λs
  results are propped up by pairs that cannot physically be generated
  (e.g. a 433 nm signal would need a 6370 nm idler).
- **`idler_theta_window=(lo, hi)`** — restricts the ξi quadrature to idlers
  whose *external* polar angle lies in the window (Snell: κ is conserved
  across the exit face). This turns `Ns` from the singles *marginal* into the
  **joint signal-AND-idler rate** — the coincidence kernel used by
  `coincidence_counts.py`. Default `None` reproduces the singles marginal
  bit-identically.
- **Per-call parameter overrides** (`lam_p, Pp, L, W, deff`): the module-level
  defaults are the paper's values (405 nm, 80 mW, 3 mm, W = √(60×30) µm ≈
  42.4 µm, 1.75 pm/V); the setup scripts override them per call, so
  paper-reproduction and setup-prediction never contaminate each other.

#### Sellmeier indices (two selectable sets)

Signal/idler use the ordinary index `n_o(λ)`; the pump uses the
angle-dependent extraordinary `n_eff(λ, θm)` (paper Eq. 4). Two coefficient
sets are available via `use_sellmeier(...)`:

- **ref 1 "kato" (default)** — Kato 1986. More accurate: agrees with the
  modern measurement (Tamošauskas 2018) to ~1–2e−4 at 405/810 nm.
- **ref 2 "eimerl"** — Eimerl 1987, what the *paper* used (via the NIST
  phase-matching program). Select it only to reproduce the paper's
  figures/labels exactly.

The difference matters near degeneracy: collinear-degenerate phase matching
sits at θm ≈ 28.82° with Kato vs 28.67° with Eimerl, so the paper's θm labels
carry Eimerl's ~0.18° offset. Every script prints the active set into its
teed log.

#### Paper errata (established in this project)

The printed **Eq. (9) is correct** as an external-angle, per-unit-azimuth
density. The genuine misprints are in the main-text **Eq. (7)**: its
prefactor should be `8π⁴` (as in the Appendix), not `2π⁴`, and its exponent
should read `exp(−½ξ²)`, not `+`. The `8π⁴` form reproduces Table I's
simulated η ≈ 2.6e−10, the independent closed form Eq. (8) (2.63e−10), and
the measured 2.7–2.8e−10. Full derivation in
`src/archive/spdc_eq9_note.tex`.

### `src/plot_fig6.py` — recreation of the paper's Fig. 6(e)–(h)

Four stacked log-scale color maps of Ns(λs, θs) — λs = 430–1000 nm,
θs = −6°…+6° — one per phase-matching angle θm = 28.6°, 28.8°, 29.1°, 29.4°,
with the perfect-phase-matching tuning curve overlaid as a white dashed line,
as in the paper. Runs the sanity check first, then writes two variants:

- `output/spdc_fig6.png` — **Eimerl** indices (the paper's set): matches the
  published panels exactly, e.g. the 28.6° panel's collinear pair at
  ~754/875 nm.
- `output/spdc_fig6_kato.png` — **Kato** indices (the library default, more
  accurate): same θm labels give visibly different maps near degeneracy.

Each map is built column-by-column (one λ at a time) to keep the memory of
the ξi quadrature axis bounded.

### `src/detector_counts.py` — singles rate on a real detector (my setup)

Absolute photons/s on a finite-aperture detector for **this lab's setup**
(not the paper's). Every knob is in the `PARAMETERS` block at the top:

- pump 405 nm, **40 mW** laser output, collimated, 2 mm beam diameter
  (W = 1 mm); BBO **L = 0.1 mm**, deff = 1.75 pm/V;
- detector at θs = **3°** (lab angle), **625 mm** from the crystal, **6 mm**
  circular opening;
- Gaussian bandpass **10 nm FWHM** centred on the degenerate 810 nm.

What it computes:

1. **Ideal cut angle** (`ideal_cut_deg`): the internal phase-matching angle
   θm whose degenerate cone passes exactly through the detector at normal
   incidence — θm ≈ 29.24° for 3°. Solved by bisection on the tuning curve.
2. **Circular-aperture geometry, exact**: `Ns` gives counts/s into the full
   2π azimuthal ring per (dλ × dθ) bin. A circular aperture (radius r,
   distance D) centred on the ring at R0 = D·tan θs intersects the cone of
   polar angle θ (ring radius ρ = D·tan θ) over the azimuthal arc
   `dφ(θ) = 2·arccos[(ρ² + R0² − r²)/(2ρR0)]`, so the detected rate is
   `∫dθ [dN/dθ]·dφ(θ)/2π` with `[dN/dθ]` the filter-weighted,
   wavelength-integrated ring density.
3. **Component efficiencies** (`EFFICIENCIES` block, datasheet values):
   the pump chain (2× F01 mirror 0.92, HWP 0.999, AR-coated crystal face)
   gives η_pump = 0.846, so 33.8 mW reach the crystal — `PP` is the *laser
   output* power and `Ns` is fed `Pp = PP·η_pump`. The collection chain
   (crystal exit face 1.0 × P01 mirror 0.958 × fiber coupling 0.80 (all-in) ×
   fiber 1.0 × APD 0.30, × filter peak 0.983) gives η_collect = 0.226,
   multiplying the detected rate. The Gaussian filter shape stays
   peak-normalised; the 0.983 peak transmission is carried separately.
4. **Crystal-rotation scan**: the detected rate vs external rotation angle α,
   with a secondary top axis in θm. The α↔θm mapping accounts for pump
   refraction at the entrance face: `sin α = n_eff(λp, θm)·sin(θm − θ_cut)` —
   a given Δθm needs a ~1.66× larger external rotation.

Headline numbers (with efficiencies): **1.24e3 /s detected** on one detector,
full azimuthal ring 5.40e4 /s (loss-free: 6.46e3 and 2.83e5 /s). The scan is a
broad, symmetric thin-crystal sinc² pattern, 2 orders down at α ≈ ±3.1°.
Writes `output/detector_scan_mysetup.png`.

### `src/coincidence_counts.py` — signal–idler coincidences (my setup)

Coincidence rate between **two** detectors at +3° and −3°, reusing
`detector_counts.py`'s pump/crystal parameters, ideal-cut solver, and
efficiency block (both arms have identical chains). Per-arm parameters
(angle, distance, aperture, filter; `filter_fwhm=None` = open arm) in the
`PARAMETERS` block.

The key idea: `Ns` is a *marginal* — its ξi quadrature integrates the idler
over everything — but its integrand `exp(−½(ξs−ξi)²)·sinc²(½LΔkz)` is the
**joint signal–idler transverse-momentum density**. The coincidence rate is
therefore the same integral with three modifications:

- **(a) polar:** the ξi quadrature restricted to arm B's polar acceptance
  (`Ns(..., idler_theta_window=...)`);
- **(b) spectral:** the wavelength grid weighted by *both* filters,
  `T_A(λs)·T_B(λi(λs))` with λi the energy-conserving conjugate;
- **(c) azimuthal:** the singles' arc factor dφ_A/2π replaced by the overlap
  of arm A's arc with arm B's arc mirrored through the pump axis, smeared by
  the azimuthal momentum-correlation width σφ = (1/W)/(k0·sin θ) ≈ 0.14°
  (tiny compared to the ~10.5° arc, so geometry costs only a percent-level
  rim loss).

Physical conclusion: for this geometry the pair correlations (polar blur
≈ 0.007°, azimuthal 0.14°, λ-walk-off ≤ 0.04°) are all far smaller than the
aperture, so the geometry is near-lossless and **the filters set the
heralding ratio**: with identical 10 nm Gaussians in both arms the spectral
factor is T², giving Rc/Rs → 1/√2 ≈ 0.71 (model: 0.689 including geometry).
Component losses multiply on top — singles by η_arm, coincidences by
η_A·η_B — so the observable heralding picks up a factor η_B:

- singles **1.24e3 /s** per arm, coincidences **192 /s**, heralding
  **0.156** = η_B·0.689 (loss-free: 6.46e3, 4.45e3, 0.689);
- accidentals R_A·R_B·τ ≈ 1.5e−3 /s at τ = 1 ns;
- built-in validation printed on every run: opening arm B completely must
  recover arm A's singles rate exactly (rel. diff 0).

Writes `output/coincidence_scan.png` (singles + coincidences + heralding vs
α). The heralding ratio is flat ≈ 0.69 across the whole scan — near-degenerate
filters keep the twins mirrored, so only the filters set Rc/Rs, not geometry.

### Shared helpers

- **`src/detector_scan_common.py`** — scan machinery shared by the detector
  and coincidence scripts (and the archived scans): `LAM_I_MAX` (3500 nm
  idler cutoff), `spectra_vs_thetam`, `rotation_axes` (α axis + θm
  secondary-axis interpolators), `two_order_crossings`, `summarize`,
  `draw_panel`. No physics of its own.
- **`src/textlog.py`** — `tee_stdout(filename)` context manager; every
  runnable script wraps its `__main__` body in it so the printed summary is
  also written to a `.txt` in `output/`.

---

## Known limitation / planned improvement: MMF collection model

The coincidence model treats each collection arm as a bare far-field angular
bucket (6 mm aperture at 625 mm) with an implicitly plane-wave pump. In that
picture the singles cannot depend on the pump spot size, and focusing the
pump can only *reduce* coincidences (the azimuthal correlation width
σφ ∝ 1/W grows and spills past the aperture; the model predicts heralding
0.156 → 0.129 for W = 1 → 0.1 mm).

**The experiment shows the opposite sign**: focusing the pump raised singles
~15%, coincidences ~32%, and heralding ~14% — the fingerprint of improved
pump–collection *spatial-mode overlap*, which the momentum-only model
structurally cannot produce. The real collection is a multimode fiber
(50 µm core, 0.22 NA) behind an f = 11 mm collimator: a finite-étendue
spatial+angular mode with a ~2.4 mm footprint at the crystal (≈ matched to
the 2 mm pump) and ±0.13° angular acceptance.

The planned fix is to replace the hard-aperture kernel with a
**finite-étendue MMF collection mode** per arm (Gaussian pump profile +
mode-overlap integrals), validate against the collimated/focused data points,
then scan the pump waist — the model predicts an optimum near ~0.1 mm pump
diameter, where pump divergence ≈ fiber acceptance. Full analysis, measured
rates, and implementation plan in **`mmf_focusing_findings.md`**.

---

## Archive

`src/archive/` holds material that is superseded or kept for reference only
(the scripts still run from there; their outputs go to `output/archive/`):

- **`plot_detector_scan.py`** / **`plot_detector_scan_bandpass.py`** — the
  original fixed-detector (3°) / rotating-crystal scans for the *paper's*
  parameters, used to establish the broadband-vs-filtered asymmetry: the
  broadband curve is strongly asymmetric (sharp degenerate-caustic edge at
  +0.8°, slow −7° tail from the tuning-curve arms), while a 10 nm bandpass
  makes it symmetric (±0.3°). Superseded by `detector_counts.py`, which
  models the real aperture and efficiencies.
- **`spdc_eq9_note.tex`** — LaTeX write-up of the paper's angle conventions
  and errata (why Eq. (9) is correct as printed and the misprints are in
  Eq. (7)). Built PDF in `output/archive/spdc_eq9_note.pdf`.
- **`code.txt`** — the original Mathematica implementation of Eq. (9) that
  this project debugged and rewrote in Python; reference only.
