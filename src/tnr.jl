# Thermal runaway integrator: couples envelope build-up, compressional
# heating, degeneracy-suppressed heat capacity, CNO burning, and convective
# dredge-up mixing into one stiff ODE, integrated continuously from the start
# of accretion through ignition to the point degeneracy lifts (the runaway
# "peak"), at which point control passes to eruption.jl for the expansion/
# mass-ejection phase (a physically distinct regime this envelope model does
# not attempt to resolve).
#
# State vector u = [M_env, T, X_H, X_CNO]:
#   M_env  : accreted envelope mass [g]
#   T      : envelope base temperature [K]
#   X_H    : base hydrogen mass fraction (depleted by burning)
#   X_CNO  : base CNO-group mass fraction (raised by convective dredge-up)
# X_He is not tracked explicitly; it is recovered as 1 - X_H - X_CNO wherever
# needed (three-species composition, consistent with envelope.jl's EOS/heat
# capacity functions).

const Q_H_TO_HE = 6.448e18  # erg/g released converting H to He (26.73 MeV per 4 amu)

struct TNRParams
    M_wd::Float64
    R_wd::Float64
    wd_composition::WDComposition
    mixing_fraction::Float64
    convection_alpha::Float64
    t_dyn::Float64
    mdot_profile::AnyMdotProfile
end

function TNRParams(p::NovaParams)
    R_wd = wd_radius(p.M_wd, p.wd_composition)
    t_dyn = sqrt(R_wd^3 / (G * p.M_wd))
    return TNRParams(p.M_wd, R_wd, p.wd_composition, p.mixing_fraction, p.convection_alpha, t_dyn, p.mdot_profile)
end

"""
    envelope_state(u, tp::TNRParams) -> (rho, mue, T_F, degeneracy_ratio)

Derived local state (density, mean molecular weight per electron, Fermi
temperature, and T/T_F) at the envelope base for ODE state `u`.
"""
function envelope_state(u, tp::TNRParams)
    M_env, T, X_H, X_CNO = u
    X_He = max(1 - X_H - X_CNO, 0.0)
    mue = mu_e(X_H, X_He, X_CNO)
    P_base = envelope_base_pressure(M_env, tp.M_wd, tp.R_wd)
    rho = envelope_base_density(P_base, mue)
    T_F = fermi_temperature(rho, mue)
    return rho, mue, T_F, T / T_F, X_He
end

smoothstep(x) = 1 / (1 + exp(-x))

function tnr_rhs!(du, u, tp::TNRParams, t)
    M_env, T, X_H, X_CNO = u
    rho, _mue, _T_F, degen_ratio, X_He = envelope_state(u, tp)

    mdot = tp.mdot_profile(t)
    du[1] = mdot

    eps_nuc = epsilon_nuc(T, rho, X_H, X_CNO)
    eps_loss = epsilon_nu(T)
    cp = specific_heat_capacity(rho, T, X_H, X_He, X_CNO)

    # Compressional (PdV) heating from the growing weight of the envelope:
    # specific heating rate (P/ρ)·d(ln ρ)/dt, with d(ln ρ)/dt = (3/5)d(ln P)/dt
    # from the degenerate EOS (P∝ρ^(5/3)) and d(ln P)/dt = ṁ/M_env (since
    # P_base ∝ M_env at fixed M_wd,R_wd). Routed through the same
    # degeneracy-suppressed c_p as nuclear heating — the smaller c_p becomes
    # as degeneracy sets in, the more effectively a fixed compressional power
    # input raises T, exactly the feedback that drives a nova.
    P_base = envelope_base_pressure(M_env, tp.M_wd, tp.R_wd)
    eps_compr = M_env > 0 ? (P_base / rho) * (3 / 5) * (mdot / M_env) : 0.0

    du[2] = (eps_nuc - eps_loss + eps_compr) / cp

    du[3] = -eps_nuc / Q_H_TO_HE

    # Convective dredge-up: activates only once the layer is closing in on
    # degeneracy lifting (T/T_F beyond ~0.85, a proxy for the onset of a
    # superadiabatic, convectively unstable gradient during the runaway
    # itself — not during quiescent accretion, where the envelope is safely
    # degenerate with T/T_F well below this), mixing in core material on an
    # order-unity multiple of the dynamical time. The transition must be
    # narrow: this factor multiplies a *dynamical* (second-scale) rate, so
    # even a seemingly negligible ~1e-6 leak, sustained over the years-long
    # accretion phase, is enough to fully saturate the mixing well before
    # the real runaway.
    conv_onset = smoothstep((degen_ratio - 0.85) / 0.01)
    X_CNO_target = (1 - tp.mixing_fraction) * SOLAR_ENVELOPE.X_CNO +
                    tp.mixing_fraction * core_cno_reservoir(tp.wd_composition)
    du[4] = tp.convection_alpha * conv_onset * (X_CNO_target - X_CNO) / tp.t_dyn

    return nothing
end

"""
    run_tnr(params::NovaParams; M_env0=1.0e-5Msun, tmax=1.0e6yr_s) -> sol, tp

Integrates the coupled accretion/thermal-runaway ODE starting from a thin but
already electron-degenerate seed envelope (`M_env0` is chosen so the Fermi
temperature at the resulting base density sits comfortably above
`core_temperature` — starting non-degenerate would make the base state
`T > T_F` from t=0, and the degeneracy-lifting termination callback below,
which looks for a `T/T_F` crossing, would never fire), through compressional
heating and CNO ignition, to the point
degeneracy lifts (T crosses T_F at the envelope base — the runaway "peak"),
where integration terminates via a `ContinuousCallback`. Uses `Rodas5`, a
stiff Rosenbrock solver, since the right-hand side spans years of slow
accretion and seconds-minutes of explosive runaway within one integration.
"""
function run_tnr(params::NovaParams; M_env0 = 1.0e-5 * Msun, tmax = 1.0e6 * yr_s)
    tp = TNRParams(params)
    u0 = [M_env0, params.core_temperature, SOLAR_ENVELOPE.X_H, SOLAR_ENVELOPE.X_CNO]

    condition(u, _t, _integrator) = envelope_state(u, tp)[4] - 1.0  # degen_ratio - 1
    affect!(integrator) = SciMLBase.terminate!(integrator)
    cb = OrdinaryDiffEq.ContinuousCallback(condition, affect!)

    prob = OrdinaryDiffEq.ODEProblem(tnr_rhs!, u0, (0.0, tmax), tp)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEqRosenbrock.Rodas5(); callback = cb,
                                 reltol = 1.0e-8, abstol = 1.0e-12 .* abs.(u0) .+ 1.0e-30,
                                 maxiters = 1_000_000)
    return sol, tp
end
