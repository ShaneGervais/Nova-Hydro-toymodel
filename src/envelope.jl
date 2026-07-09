# Envelope equation of state, hydrostatic base state, and nuclear energy
# generation for the accreted layer on the WD surface.
#
# Simplification note: physics here follows a *single representative shell*
# at the envelope base (the classic approach of Fujimoto 1982's analytic
# runaway analysis) rather than a fully resolved multi-zone radiative-diffusion
# stellar structure. `NovaParams.n_zones` is used only by `visualization.jl` to
# interpolate a plausible radial T/ρ profile shape for the movie's envelope
# cutaway, not to resolve independent physics per shell. See README for the
# full list of documented simplifications.

const SOLAR_ENVELOPE = (X_H = 0.70, X_He = 0.28, X_CNO = 0.02)

"""
    core_cno_reservoir(wd_composition) -> Float64

Fraction of dredged-up WD core material that acts as CNO-cycle catalyst once
mixed into the H-rich envelope. CO-core material (carbon+oxygen) is itself
CNO-cycle fuel, so essentially all of it counts; ONe-core material includes
some Ne/Mg that instead feeds the (slower, less energetic) Ne-Na/Mg-Al chain,
approximated here as a modest reduction in effective catalytic yield.
"""
core_cno_reservoir(wd_composition::WDComposition) = wd_composition == CO ? 1.0 : 0.85

"""
    mu_e(X_H, X_He, X_CNO) -> mean molecular weight per electron

Assumes full ionization: 1/μ_e = X_H·1 + X_He·(1/2) + X_CNO·(1/2) (CNO-group
nuclei all have Z/A ≈ 1/2, like He).
"""
mu_e(X_H, X_He, X_CNO) = 1 / (X_H * 1.0 + X_He * 0.5 + X_CNO * 0.5)

"""
    mu_ion(X_H, X_He, X_CNO) -> mean molecular weight per ion

1/μ_ion = X_H + X_He/4 + X_CNO/14 (CNO group approximated as A≈14).
"""
mu_ion(X_H, X_He, X_CNO) = 1 / (X_H + X_He / 4 + X_CNO / 14)

"""
    mu_total(X_H, X_He, X_CNO) -> mean molecular weight (ions + electrons)

Standard fully-ionized-gas mean molecular weight: 1/μ = 2X + 3Y/4 + Z/2.
"""
mu_total(X_H, X_He, X_CNO) = 1 / (2 * X_H + 0.75 * X_He + 0.5 * X_CNO)

# --- Equation of state --------------------------------------------------

# Non-relativistic degenerate electron pressure coefficient, built from
# fundamental constants: P_deg = K_DEG_NR * (ρ/μ_e)^(5/3).
const K_DEG_NR = hbar^2 * (3 * pi^2)^(2 / 3) / (5 * m_e * m_u^(5 / 3))

"""
    degenerate_pressure(rho, mue) -> P [erg/cm^3]

Non-relativistic zero-temperature degenerate electron pressure. Valid for the
moderately dense, non-relativistic WD envelope base considered here (breaks
down only near the Chandrasekhar mass, out of scope for classical novae).
"""
degenerate_pressure(rho, mue) = K_DEG_NR * (rho / mue)^(5 / 3)

"""
    fermi_temperature(rho, mue) -> T_F [K]

Electron Fermi temperature T_F = E_F/k_B with the non-relativistic Fermi
energy E_F = (ħ²/2mₑ)(3π² n_e)^(2/3). Sets the scale below which electrons are
degenerate and their heat capacity is suppressed — the actual thermal-runaway
mechanism in a nova (Kippenhahn & Weigert, *Stellar Structure and Evolution*,
ch. 15).
"""
function fermi_temperature(rho, mue)
    n_e = rho / (mue * m_u)
    E_F = hbar^2 / (2 * m_e) * (3 * pi^2 * n_e)^(2 / 3)
    return E_F / k_B
end

ideal_gas_pressure(rho, T, mu) = rho * k_B * T / (mu * m_u)
radiation_pressure(T) = a_rad * T^4 / 3

"""
    total_pressure(rho, T, X_H, X_He, X_CNO) -> P [erg/cm^3]

Sum of degenerate electron, ideal-gas, and radiation pressure — a documented
simplification of the true partially-degenerate equation of state (a real
stellar-structure code solves the full Fermi-Dirac integrals; summing limits
is accurate to within a factor of order unity in the partially-degenerate
transition region, which is adequate for driving this simulation's dynamics).
"""
function total_pressure(rho, T, X_H, X_He, X_CNO)
    mue = mu_e(X_H, X_He, X_CNO)
    mu = mu_total(X_H, X_He, X_CNO)
    return degenerate_pressure(rho, mue) + ideal_gas_pressure(rho, T, mu) + radiation_pressure(T)
end

"""
    specific_heat_capacity(rho, T, X_H, X_He, X_CNO) -> c_p [erg g^-1 K^-1]

Ions always contribute classically (they are non-degenerate at these
densities); electrons contribute their classical 3/2 k_B/(μ_e m_u) only once
T approaches the Fermi temperature, tapering to ~0 as T ≪ T_F. This
degeneracy-suppressed heat capacity is *why* a nova runs away: heating a
degenerate layer barely raises its pressure (no compensating expansion/
cooling) because most of the added energy has nowhere to go until degeneracy
lifts.
"""
function specific_heat_capacity(rho, T, X_H, X_He, X_CNO)
    mue = mu_e(X_H, X_He, X_CNO)
    mui = mu_ion(X_H, X_He, X_CNO)
    T_F = fermi_temperature(rho, mue)
    c_ion = 1.5 * k_B / (mui * m_u)
    c_e = 1.5 * k_B / (mue * m_u) * min(1.0, T / T_F)
    return c_ion + c_e
end

# --- Hydrostatic envelope base state -------------------------------------

"""
    envelope_base_pressure(M_env, M_wd, R_wd) -> P_base [erg/cm^3]

Thin-envelope hydrostatic estimate: P_base ≈ g·Σ with surface gravity
g = G M_wd/R_wd² and column density Σ = M_env/(4π R_wd²), i.e.
P_base = G M_wd M_env / (4π R_wd⁴). Standard in nova envelope models (e.g.
Fujimoto 1982; Prialnik & Kovetz 1995) where the accreted layer is thin
compared to the WD radius.
"""
envelope_base_pressure(M_env, M_wd, R_wd) = G * M_wd * M_env / (4 * pi * R_wd^4)

"""
    envelope_base_density(P_base, mue) -> rho [g/cm^3]

Inverts the non-relativistic degenerate EOS for density, appropriate while
degenerate pressure dominates the envelope base (true for essentially the
whole accretion phase up to ignition).
"""
envelope_base_density(P_base, mue) = mue * (P_base / K_DEG_NR)^(3 / 5)

# --- Nuclear energy generation --------------------------------------------

const T_HCNO_ONSET = 8.0e7  # K, approximate cold->hot CNO transition (Iliadis, "Nuclear Physics of Stars")

"""
    epsilon_nuc(T, rho, X_H, X_CNO) -> erg g^-1 s^-1

CNO-cycle hydrogen burning rate, blended between two physically distinct
regimes:
- Cold/normal CNO (T ≲ 8×10⁷ K): strongly temperature-sensitive,
  capture-rate-limited, ε ∝ ρ X_H X_CNO T^18 (illustrative power-law index in
  the literature range 13-20 for this temperature window).
- Hot CNO (T ≳ 8×10⁷ K): β-decay-limited (rate-limiting steps ¹⁴O, ¹⁵O with
  half-lives ~70-120 s), essentially independent of ρ and T:
  ε ≈ (Q/τ_β)·(X_CNO/A m_u) ~ 8×10¹⁵·X_CNO erg g⁻¹ s⁻¹, derived here from
  Q≈25 MeV per catalytic cycle and τ_β≈200 s.

The two are joined with a smooth logistic in T so the ODE integrator sees a
continuous right-hand side; normalizations are order-of-magnitude/continuity
calibrated, not a precision nuclear-reaction-network result (see README).
"""
function epsilon_nuc(T, rho, X_H, X_CNO)
    T7 = T / 1.0e7
    eps_cold = 2.0e-4 * rho * X_H * X_CNO * T7^18
    eps_hot = 8.0e15 * X_CNO
    # Narrow transition width: eps_hot is ~15-17 orders of magnitude larger
    # than eps_cold far below the transition, so even a seemingly negligible
    # blend weight (e.g. 1e-4) would leak an unphysically large hot-CNO
    # contribution at temperatures where that channel cannot actually operate.
    # A width of a few % of T_HCNO_ONSET keeps blend at true machine-zero
    # there while still giving the ODE solver a smooth crossing.
    blend = 1 / (1 + exp(-(T - T_HCNO_ONSET) / (0.02 * T_HCNO_ONSET)))
    return (1 - blend) * eps_cold + blend * eps_hot
end

"""
    epsilon_nu(T) -> erg g^-1 s^-1

Illustrative order-of-magnitude thermal-neutrino loss term, negligible below
~2×10⁸ K and only mildly relevant near peak temperature; not a precision
Itoh et al. (1996) fit.
"""
epsilon_nu(T) = 1.0e-2 * (T / 3.0e8)^9
