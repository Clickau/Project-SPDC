# SPDC-in-BBO flux calculations

Numerical recreation of the angle-resolved parametric-fluorescence flux from
**Hsu & Lai, "Angle-Resolved Spectroscopy of Parametric Fluorescence"
(arXiv:1302.6226)**. Type-I (e → o + o) SPDC in BBO, pump λp = 405 nm,
degenerate signal/idler at 810 nm. The user is a physicist; treat the physics
as the source of truth and keep the scripts numerically consistent with each
other (they share `spdc_physics.py` — **edit physics there, nowhere else**).

`paper.pdf` is the reference paper (equations verified against the rendered PDF,
not just text extraction). `code.txt` is the user's original Mathematica
implementation of Eq.(9) that was debugged/rewritten in Python (see bug list
below); keep it only for reference.

## Files (post 2026-07-04 refactor: physics library + thin plot scripts)

- **`spdc_physics.py`** — the physics library; all constants and physics live
  here, the plot scripts contain none. Defines the SI constants, experiment
  parameters, BBO Sellmeier (`n_o`, `n_e_principal`, `n_eff`), `sinc`, `Np`,
  the unified flux `Ns(lam_s, theta_s, ...)` (numpy-broadcasts λs against θs,
  so one function serves both the Fig. 6 maps and the detector scans; keyword
  options: `external` angle convention, `Nxi` quadrature points, `lam_i_max`
  idler cutoff), `tuning_curve`, and `crystal_rotation`. Running
  `python spdc_physics.py` executes `sanity_check()` against the paper.
- **`detector_scan_common.py`** — shared machinery for the two detector-scan
  diagrams: `LAM_I_MAX` (3500 nm idler cutoff), `spectra_vs_thetam` (the
  θm × λ count matrix every panel derives from), `rotation_axes` (α axis +
  θm secondary-axis interpolators), `two_order_crossings`, `summarize`,
  `draw_panel`.
- **`plot_fig6.py`** — Fig. 6 recreation. Writes `spdc_fig6.png`: the four
  stacked log-scale λs–θs maps Fig. 6(e)–(h) with white-dashed
  perfect-phase-matching tuning curves. Also runs the sanity check first.
- **`plot_detector_scan.py`** — fixed detector at external angle θs = 3°,
  rotate the crystal, plot wavelength-integrated count vs crystal rotation
  angle α. Two panels: broadband + hard 805–815 nm top-hat. Writes
  `detector_scan.png`.
- **`plot_detector_scan_bandpass.py`** — same scan with reworked wavelength
  selection, three panels: (1) broadband (all λ), (2) Gaussian bandpass FWHM
  10 nm @ 810 nm, (3) monochromatic spectral density dN/dλ at 810 nm
  (counts/s/nm). Writes `detector_scan_bandpass.png`.
- **`spdc_eq9_note.tex`** (`.log`) — derivation and write-up of the two Eq.(9)
  typo corrections (see below).

(History: the pre-refactor layout had the physics inside `spdc_fig6.py` and a
duplicated flux kernel in each detector script. The unified `Ns` was verified
**bitwise identical** to both old kernels — and `tuning_curve` /
`crystal_rotation` / the sanity totals likewise — before the old files were
deleted.)

## Master equation — Eq.(9), angle-resolved signal flux

As implemented in `spdc_physics.py` (the **corrected** normalisation):

```
Ns(λs,θs) = 2π·√(2π) · ħ·deff²·ωp·L²·Np / (ε0·ns·ni·np·λs⁵·λi)
            · sin(2θs) · dλs · dθs
            · ∫dξi  exp(−½(ξs−ξi)²) · sinc²(½·L·Δkz)
```

- `sinc(x) = sin(x)/x`, so `sinc²` is squared (NumPy's `np.sinc` is
  `sin(πx)/(πx)` — divide the argument by π).
- ξs = ks·sin(θs)·W, ξi = κi·W (dimensionless transverse momenta);
  W = √(60 µm × 30 µm) ≈ 42.4 µm.
- Longitudinal mismatch (Appendix form):
  `Δkz = ks·cos(θs) + √(ki² − (ξi/W)²) − kp`, with `k = 2π·n/λ` for each wave.
- Energy conservation → idler: `λi = λs·λp/(λs − λp)`.

### Two independent typos in the paper's printed Eq.(9) — IMPORTANT

The printed Eq.(9) is **~2π·ns² ≈ 17× too low** vs the paper's own Table I /
Fig. 6 (gives η ≈ 1.6e-11 instead of ≈ 2.7e-10). The corrected form above is
what reproduces the paper. Both typos are derived in `spdc_eq9_note.tex`:

1. **Factor ~17 (= 2π·ns²):** dropped in the appendix→Eq.(9) substitution
   κs = ks·sin θs (the κs·dκs angular-Jacobian step). Correct Jacobian
   `κs·dκs = ½·ks²·sin2θs·dθs` with `ks = 2π·ns/λs`. The printed Eq.(9) keeps
   the 1/λs² scaling (the λs⁵ in the denominator) but drops the numeric+index
   coefficient `2π·ns²`. The ns² fingerprint can only come from the in-medium
   wavenumber ks, so the slip is entirely in this step.
2. **Factor 4:** Eq.(7)'s prefactor is `2π⁴` in the main text but `8π⁴` in the
   Appendix. `8π⁴` is correct (gives η = 2.72e-10); `2π⁴` is the misprint. The
   relabel ∫d²ξ→∫d²ξi is a unit-Jacobian shift and cannot produce the 4.
   (Minor third typo: Eq.(7) exponent prints `exp(+½ξ²)`; must be `−½ξ²`, as the
   Appendix has it.)

## Indices, Sellmeier, constants

- Signal & idler ordinary: `ns = no(λs)`, `ni = no(λi)`. Pump effective
  extraordinary (Eq. 4): `np = n_eff(θm,λp) = [cos²θm/no² + sin²θm/ne²]^(−½)`.
- BBO Sellmeier (λ in µm):
  - `no² = 2.7359 + 0.01878/(λ²−0.01822) − 0.01354·λ²`
  - `ne² = 2.3753 + 0.01224/(λ²−0.01667) − 0.01516·λ²`
- Constants/params: ħ = 1.054571817e−34 J·s; ε0 = 8.8541878188e−12 F/m;
  c = 2.99792458e8 m/s; deff = 1.75 pm/V; L = 3 mm; λp = 405 nm; Pp = 80 mW ⟹
  Np = Pp/(ħ·ωp) = 1.63e17 /s.
- **Angle convention:** Eq.(9) internally uses the *internal* angle
  (κs = ks·sinθ, ks = ns·2π/λs); Fig. 6 plots the *external* (air) angle
  (κs = k0·sinθ_ext, k0 = 2π/λs). Both give identical total flux — only the
  y-axis scale differs by ~ns. Use external for figures.

## Sanity checks (must hold)

- η = 2Ns/Np ≈ 2.6–2.7e−10 at 810 nm for θm ≈ 28.8–29.4°.
- Total degenerate flux 2Ns ≈ 4.2e7 /s at Pp = 80 mW (code gives 4.44e7).
- Fig. 6 phase-matching angles θm = 28.6°, 28.8°, 29.1°, 29.4°;
  λs ≈ 430–1000 nm, θs ≈ −6°…+6°, color on log scale.
- Note: collinear-degenerate (810 nm at θ=0) lands at θm ≈ 28.82° with this
  Sellmeier vs the paper's 28.6° (paper used a NIST index program) — a ~0.2°
  label offset; the physics is correct.

## Detector-scan analysis — key physical findings

Detector fixed at external θs = 3°, crystal rotated. Because the pump refracts
at the entrance face, the external crystal rotation α ≠ internal Δθm:
`sin α = n_eff(λp,θm)·sin(θm − θ_cut)` — a given Δθm needs ~n_eff ≈ 1.66× larger
external α. Plot x-axis = α, secondary top axis = θm.

- **Cone-axis tilt cancels (verified by 2D ray-trace):** rotating the crystal
  tilts the internal pump, but pump-in + signal-out refraction gives an internal
  signal–pump angle ≈ 1.80° for *all* α (1.79–1.81° over α = +0.4…−8.3°),
  because n_o(signal) ≈ n_eff(pump) ≈ 1.66 in BBO. So the simple
  κs = k0·sin(3°) model is geometrically correct.
- **Idler cap:** wavelength integration zeroes signals whose idler exceeds
  ~3500 nm (BBO transparency / index-data limit). Without it the short-λ tail is
  propped up by pairs that cannot be generated (e.g. 433 nm → 6370 nm idler).
- **Broadband asymmetry is real (broadband vs narrowband):** broadband peak at
  α ≈ +0.37° (θm ≈ 29.22°); dropping 2 orders takes only +0.8° on the high side
  (sharp degenerate-caustic edge) but −7° on the low side (the tuning-curve arms
  keep crossing 3°). A 10 nm bandpass @810 nm makes it **symmetric** (±0.3° for
  2 orders) — that is the "intuitive"/filtered case (paper Fig. 3).
- The detector-acceptance factor (dθ, azimuth) only rescales the y-axis; it does
  not change the shape or the orders-of-magnitude conclusions.
