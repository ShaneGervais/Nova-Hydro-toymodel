# Post-runaway eruption: once degeneracy lifts, the envelope becomes
# radiation-pressure-dominated and expands. Mass loss is treated with the
# Kato & Hachisu (1994) "optically thick wind" picture: photospheric
# luminosity self-regulates near the Eddington limit while the excess energy
# above the envelope's binding energy drives an outflow.

"""
    electron_scattering_opacity(X_H) -> kappa_es [cm^2/g]

κ_es = 0.2(1+X_H), the standard fully-ionized Thomson-scattering opacity.
"""
electron_scattering_opacity(X_H) = 0.2 * (1 + X_H)

"""
    eddington_luminosity(M_wd, X_H) -> L_Edd [erg/s]
"""
eddington_luminosity(M_wd, X_H) = 4 * pi * G * M_wd * c_light / electron_scattering_opacity(X_H)

envelope_binding_energy(M_wd, M_env, R_wd) = G * M_wd * M_env / R_wd

"""
    envelope_thermal_energy(T, M_env, rho) -> E [erg]

Radiation energy content of the envelope shell once non-degenerate,
E_rad = (a T^4/3) · V with shell volume V = M_env/ρ.
"""
envelope_thermal_energy(T, M_env, rho) = a_rad * T^4 / 3 * M_env / rho

const EXPLOSIVE_BURN_FRACTION = 0.1  # fraction of remaining envelope H burned during the brief post-degeneracy runaway peak

"""
    ejecta_velocity(T_peak, M_env, rho_peak, M_wd, R_wd, X_H_peak) -> v_ej [cm/s]

Energy balance: available energy in excess of the envelope's own
gravitational binding energy becomes ejecta kinetic energy. The dominant
term is *not* the radiation energy already stored at the instant degeneracy
lifts (that alone is typically far too small once one-zone modeling limits
how hot the envelope gets — see README) but the nuclear energy released as
the still-abundant unburned hydrogen ignites explosively during the brief,
unresolved expansion between degeneracy lifting and this eruption model
taking over. `EXPLOSIVE_BURN_FRACTION` (illustrative, not a precision yield)
is the fraction of the remaining envelope hydrogen assumed to burn in that
interval:

    v_ej = √(2 max(E_rad + f_burn·X_H·M_env·Q_H_TO_HE − E_bind, 0)/M_env)
"""
function ejecta_velocity(T_peak, M_env, rho_peak, M_wd, R_wd, X_H_peak)
    E_th = envelope_thermal_energy(T_peak, M_env, rho_peak)
    E_burst = EXPLOSIVE_BURN_FRACTION * X_H_peak * M_env * Q_H_TO_HE
    E_bind = envelope_binding_energy(M_wd, M_env, R_wd)
    return sqrt(2 * max(E_th + E_burst - E_bind, 0.0) / M_env)
end

"""
    wind_mass_loss_rate(L, v_ej, M_wd, R_wd) -> Mdot_wind [g/s]

Energy-conservation estimate: radiated power (capped at L_Edd) supplies both
the kinetic energy and the work to unbind the outflowing mass,
L = Ṁ_wind (½v_ej² + GM_wd/R_wd).
"""
wind_mass_loss_rate(L, v_ej, M_wd, R_wd) = L / (0.5 * v_ej^2 + G * M_wd / R_wd)

"""
    EruptionState

One frame of the post-runaway eruption/decline light curve.
"""
struct EruptionState
    t::Float64          # s, since eruption onset
    L::Float64           # erg/s
    R_photosphere::Float64  # cm
    v_ej::Float64            # cm/s
    M_env::Float64            # g, remaining envelope mass
end

"""
    run_eruption(M_wd, R_wd, T_peak, M_env_peak, rho_peak, X_H_peak; n=200) -> Vector{EruptionState}

Builds the eruption light curve: a plateau near `L_Edd` while the optically
thick wind ejects the envelope (duration ≈ M_env_peak/Ṁ_wind), then an
exponential decline (time constant `3×` the plateau duration — a simple,
documented stand-in for the real recession of the photosphere through the
optically thinning ejecta) as the remaining envelope mass empties out.
"""
function run_eruption(M_wd, R_wd, T_peak, M_env_peak, rho_peak, X_H_peak; n = 200)
    L_edd = eddington_luminosity(M_wd, X_H_peak)
    v_ej = ejecta_velocity(T_peak, M_env_peak, rho_peak, M_wd, R_wd, X_H_peak)
    mdot_wind = wind_mass_loss_rate(L_edd, v_ej, M_wd, R_wd)
    t_plateau = M_env_peak / mdot_wind
    tau_decline = 3 * t_plateau
    t_end = t_plateau + 6 * tau_decline

    # Log-spaced, not linear: the photosphere crosses any given visual scale
    # (e.g. the binary's orbital separation) in seconds-minutes at these
    # ejecta velocities, while the full plateau+decline lasts hours-days —
    # linear sampling would put almost every frame in the "already expanded
    # past the frame" tail and none where the expanding shell is still
    # visible on screen.
    states = EruptionState[]
    for t in exp.(range(log(1.0), log(t_end); length = n))
        if t <= t_plateau
            L = L_edd
            M_env = max(M_env_peak - mdot_wind * t, 0.0)
        else
            L = L_edd * exp(-(t - t_plateau) / tau_decline)
            M_env = max(M_env_peak - mdot_wind * t_plateau, 0.0) * exp(-(t - t_plateau) / tau_decline)
        end
        R_phot = R_wd + v_ej * t
        push!(states, EruptionState(t, L, R_phot, v_ej, M_env))
    end
    return states
end
