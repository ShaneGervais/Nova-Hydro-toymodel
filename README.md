# NovaSim.jl

A visual, physically-grounded simulation of a classical nova: Roche-lobe
overflow mass transfer in a binary → thermonuclear runaway (TNR) on a white
dwarf's surface → eruption and mass ejection, rendered as a single movie.
Built to be tinkered with — every physical quantity a nova textbook discusses
(accretion rate and its time-dependence, white dwarf mass/composition,
companion mass/type, mixing, convection, densities, orbital/Roche geometry...)
is a real parameter, not a cosmetic slider.

This is Phase 1 of a two-phase project: a 1D/semi-analytic model (this
package). Phase 2 — a full multi-D hydrodynamic solver — is future work; see
[Extending to multi-D hydrodynamics](#extending-to-multi-d-hydrodynamics)
below for how it plugs in.

## Quick start

`NovaSim` itself only depends on `Makie` (the backend-agnostic plotting API),
the ODE solver, and `Roots` — it deliberately does *not* bundle a concrete
Makie backend, so loading it never requires a display/GPU. You choose the
backend (see [Backend choice](#backend-choice-glmakie-vs-cairomakie) below)
and add it yourself, once:

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.add("GLMakie")   # or "CairoMakie" on a headless machine — your choice

using GLMakie   # or `using CairoMakie`
using NovaSim

sol = run_nova(NovaParams())   # the "common" classical nova default (see below)
make_movie(sol, "nova.mp4")
```

or just run `examples/classical_nova_co_wd.jl` (after `Pkg.add`-ing a backend).

### The default scenario

`NovaParams()` with no arguments is a typical classical nova / cataclysmic
variable: a 0.8 M_sun CO white dwarf accreting at 5×10⁻¹⁰ M_sun/yr from a
0.4 M_sun main-sequence donor filling its Roche lobe (orbital period a few
hours, set automatically by Roche geometry). Change anything:

```julia
NovaParams(
    M_wd = 1.0Msun, wd_composition = ONe,           # more massive ONe WD
    M_companion = 0.3Msun, companion_type = Subgiant,
    mdot_profile = FlickeringMdot(1e-9Msun/yr_s, 3e-10Msun/yr_s, 3600.0),
    mixing_fraction = 0.5, convection_alpha = 2.0,
)
```

## What each module does

| File | Physics |
|---|---|
| `orbit.jl` | Kepler's third law; Eggleton (1983) Roche lobe radius; L1 point and full Roche-lobe equipotential contour (root-found from the actual corotating-frame potential, not hand-drawn); WD (Nauenberg 1972) and companion mass-radius relations |
| `accretion.jl` | Circularization radius (from conservation of specific angular momentum at L1); disk-vs-direct-impact criterion; ballistic L1 stream trajectory (restricted 3-body integration); steady-state Shakura & Sunyaev (1973) α-disk |
| `envelope.jl` | Degenerate + ideal + radiation pressure EOS; degeneracy-suppressed heat capacity (the actual thermal-runaway mechanism); hydrostatic envelope base state; CNO/hot-CNO nuclear energy generation |
| `tnr.jl` | The coupled ODE: envelope build-up, compressional heating, nuclear burning, convective dredge-up mixing — integrated continuously from first accretion through ignition to the point degeneracy lifts |
| `eruption.jl` | Kato & Hachisu (1994) optically-thick-wind picture: Eddington-limited luminosity, energy-balance ejecta velocity, wind mass loss, light-curve decline |
| `simulation.jl` | Orchestration + the "movie clock" that remaps years of accretion, minutes of runaway, and days of eruption onto a fixed-fps movie |
| `visualization.jl` / `recording.jl` | The `Makie` scene and movie export (backend-agnostic — see below) |

## Documented simplifications

This is an educational tool, not a research-grade stellar-evolution code. Every
simplification below is a deliberate scope decision, not an oversight:

- **Single representative envelope shell**, not a resolved multi-zone
  radiative-diffusion structure (the classic approach of Fujimoto 1982's
  analytic runaway analysis). `NovaParams.n_zones` only controls how many
  concentric shells `visualization.jl` draws in the envelope cutaway
  (interpolated from the one physical temperature via `zone_temperatures`),
  not independent physics per zone.
- **No empirical "critical ignition mass" fit** (e.g. Yaron et al. 2005) is
  hard-coded — the ignition point emerges from directly integrating the
  coupled ODE (envelope build-up, compressional heating, degeneracy-suppressed
  heat capacity, temperature-sensitive CNO burning) until it runs away on its
  own. `find_ignition_time` in `simulation.jl` just labels where that
  happened for the movie clock, using the local e-folding timescale of T.
- **`T_peak` is the temperature at which the single-zone envelope's
  degeneracy lifts** (`tnr.jl`'s termination point), not the ~2-4×10⁸ K often
  quoted from multi-zone nova codes. A one-zone model can't resolve the
  localized hot sub-layer that lets a resolved envelope reach that much
  higher figure; for the default scenario this model reaches a few×10⁷ K
  instead — a known, documented limitation of the one-zone approximation, not
  a numerical error.
- **Ejecta velocity is powered mainly by explosive hydrogen burning**, not by
  the radiation energy already stored at the instant degeneracy lifts — at
  the one-zone `T_peak` above, that stored energy alone is far too small to
  unbind the envelope. `eruption.jl` instead assumes a fraction
  (`EXPLOSIVE_BURN_FRACTION = 0.1`, illustrative, not a precision yield) of
  the still-abundant unburned envelope hydrogen ignites during the brief,
  unresolved expansion between degeneracy lifting and the eruption model
  taking over — physically the same mechanism (Starrfield et al.) that
  actually drives real nova ejection.
- **Partially-degenerate EOS approximated as a sum of limits**
  (degenerate + ideal + radiation pressure) rather than solving the full
  Fermi-Dirac integrals a real stellar-structure code would.
- **CNO energy generation rate** uses an illustrative power law (T¹⁸) in the
  cold regime and a first-principles order-of-magnitude estimate
  (Q≈25 MeV / τ_β≈200 s) in the beta-decay-limited hot-CNO regime, blended
  smoothly at the literature-typical ~8×10⁷ K transition — not a nuclear
  reaction network. No isotope-by-isotope nucleosynthesis is tracked.
- **No orbital evolution**: separation/period are fixed for the whole movie
  (real systems slowly evolve via angular momentum loss and mass transfer;
  out of scope here).
- **Ṁ(t) is a direct user input** (`ConstantMdot`/`RampMdot`/`FlickeringMdot`,
  or any `t -> g/s` function), not derived from donor stellar evolution —
  this is intentional (it's exactly the knob requested) rather than a
  limitation.
- **Accretion disk is quasi-steady-state**, not a solved viscous-diffusion
  PDE — valid because the viscous timescale is much shorter than the
  timescale over which Ṁ(t) varies in this model.
- **Envelope cutaway thickness is exaggerated** ~1000× for visibility; real
  nova envelopes are a geometrically thin skin (~10⁻⁴ R_wd) that would be
  invisible drawn to scale.
- **Neutrino losses** (`epsilon_nu`) are an illustrative order-of-magnitude
  term, not an Itoh et al. (1996)-precision fit; they are not what terminates
  the runaway in this model (degeneracy lifting is).

## Extending to multi-D hydrodynamics (Phase 2, not implemented)

`simulation.jl` defines `AbstractNovaSolution` specifically as a pluggability
seam: anything exposing the same `.frames`/`.geometry` shape that
`visualization.jl` consumes can replace `NovaSolution` without touching the
rendering code. A future `HydroNovaSolution` — e.g. a 2D finite-volume Euler
solver (or SPH) for the mass-transfer stream/disk and the convective envelope
during runaway — would slot in there.

## Backend choice (GLMakie vs. CairoMakie)

`visualization.jl` is written against the backend-agnostic `Makie` API — it
never imports a concrete backend. Load and activate whichever one you want
before calling `make_movie`:

- `using GLMakie` — GPU-accelerated, interactive, the default for real use.
  Requires an actual display/GPU context; it will not even *load* headlessly
  (its GLFW dependency initializes a window at `using`-time, not just at
  render-time).
- `using CairoMakie` — CPU vector renderer, works headlessly (CI, servers
  without a display), used to verify this package's own plotting code in
  environments without a display.

## Testing

```julia
using Pkg; Pkg.activate("."); Pkg.test()
```

Covers the Roche lobe formula (including a dimensional cross-check via
`Unitful.jl`), Kepler mechanics, WD mass-radius relation, L1/disk geometry,
nuclear rate regime trends, and a full end-to-end smoke test of `run_nova`
against broad literature sanity ranges for a classical CO nova.

## References

- Eggleton, P. P. 1983, ApJ, 268, 368 (Roche lobe radius)
- Shakura, N. I. & Sunyaev, R. A. 1973, A&A, 24, 337 (α-disk)
- Nauenberg, M. 1972, ApJ, 175, 417 (WD mass-radius relation)
- Fujimoto, M. Y. 1982, ApJ, 257, 752 (analytic nova thermal-runaway theory)
- Kato, M. & Hachisu, I. 1994, ApJ, 437, 802 (optically thick wind)
- Iliadis, C. 2007, *Nuclear Physics of Stars* (CNO/hot-CNO cycle rates)
- Kippenhahn, R. & Weigert, A. 1990, *Stellar Structure and Evolution*
  (degenerate equation of state, heat capacity)
- Frank, J., King, A. & Raine, D. 2002, *Accretion Power in Astrophysics*
  (circularization radius, disk truncation)
- Warner, B. 1995, *Cataclysmic Variable Stars* (CV disk sizes, donor types)
