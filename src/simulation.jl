# Orchestrates a full nova cycle: static binary/Roche geometry, the coupled
# accretion+runaway ODE (tnr.jl), and the eruption light curve (eruption.jl),
# then remaps physical time onto a fixed-fps "movie clock" that gives each
# phase a fixed screen-time budget (accretion lasts years, runaway
# minutes-hours, eruption days-weeks — see plan/README).
#
# `AbstractNovaSolution` is the pluggability seam for a future multi-D hydro
# backend (Phase 2, not implemented here): anything implementing the same
# `.frames` / `.geometry` interface can be dropped into `visualization.jl`
# unchanged.

abstract type AbstractNovaSolution end

"""
    Geometry

Static binary/Roche/disk geometry for one `NovaParams` set. Orbital
separation and stellar radii are held fixed over the movie (no orbital
evolution is modeled — a documented simplification; see README) so this is
computed once per simulation rather than per frame.
"""
struct Geometry
    M_wd::Float64
    M_companion::Float64
    separation::Float64
    R_wd::Float64
    R_companion::Float64
    wd_pos::Tuple{Float64, Float64}
    companion_pos::Tuple{Float64, Float64}
    roche_xs::Vector{Float64}
    roche_ys::Vector{Float64}
    stream_xs::Vector{Float64}
    stream_ys::Vector{Float64}
    disk_forms::Bool
    disk_r_in::Float64
    disk_r_out::Float64
end

function build_geometry(params::NovaParams)
    a = isnan(params.separation) ?
        separation_for_roche_filling(params.M_wd, params.M_companion, params.companion_type) :
        params.separation
    R_wd = wd_radius(params.M_wd, params.wd_composition)
    R_comp = companion_radius(params.M_companion, params.companion_type)
    wd_pos = wd_position(params.M_wd, params.M_companion, a)
    comp_pos = companion_position(params.M_wd, params.M_companion, a)
    rxs, rys = roche_lobe_contour(params.M_wd, params.M_companion, a)
    sxs, sys, _ = stream_trajectory(params.M_wd, params.M_companion, a)
    forms = disk_forms(params.M_wd, params.M_companion, a, params.wd_composition)
    r_out = forms ? disk_outer_radius(params.M_wd, params.M_companion, a) : R_wd
    return Geometry(params.M_wd, params.M_companion, a, R_wd, R_comp, wd_pos, comp_pos,
                     rxs, rys, sxs, sys, forms, R_wd, r_out)
end

"""
    MovieFrame

One rendered instant: which physical `phase` it belongs to, the physical time
elapsed since the simulation began (`t_phys`), the compressed movie-clock time
(`t_movie`, seconds into the output video), and the physical state needed to
draw it.
"""
struct MovieFrame
    phase::Symbol            # :accretion, :runaway, or :eruption
    t_phys::Float64            # s since simulation start
    t_movie::Float64            # s into the rendered movie
    M_env::Float64                # g
    T::Float64                     # K (base T during accretion/runaway; photospheric Teff during eruption)
    rho::Float64                    # g/cm^3
    X_H::Float64
    X_CNO::Float64
    L::Float64                # erg/s
    mdot::Float64               # g/s
    R_photosphere::Float64       # cm
end

struct NovaSolution <: AbstractNovaSolution
    params::NovaParams
    geometry::Geometry
    frames::Vector{MovieFrame}
    t_ignition::Float64   # s since simulation start
    T_peak::Float64
    L_peak::Float64
    v_ejecta::Float64
    fps::Int
end

"""
    find_ignition_time(sol; efold_threshold=3600.0) -> t [s]

Locates the moment the accretion+runaway trajectory transitions from slow
compressional heating to explosive burning, defined as the first time the
local e-folding timescale of T (T / (dT/dt)) drops below `efold_threshold`
(default: 1 hour). Purely a label for the movie clock / diagnostics — the ODE
integration itself is one continuous system across both regimes.
"""
function find_ignition_time(sol; efold_threshold = 3600.0)
    t0, t1 = sol.t[1], sol.t[end]
    ts = exp.(range(log(max(t0, 1.0)), log(t1); length = 4000))
    for t in ts
        T = sol(t)[2]
        dt = t * 1.0e-6 + 1.0e-3
        dTdt = (sol(min(t + dt, t1))[2] - sol(max(t - dt, t0))[2]) / (2dt)
        if dTdt > 0 && T / dTdt < efold_threshold
            return t
        end
    end
    return t1 * 0.999  # fallback: never met threshold cleanly, use near-end
end

"""
    run_nova(params::NovaParams; fps=30, accretion_screen_s=8.0,
             runaway_screen_s=6.0, eruption_screen_s=10.0) -> NovaSolution

Runs the full accretion -> thermonuclear runaway -> eruption cycle for
`params` and packages it into per-frame data ready for `visualization.jl`,
using a fixed-fps movie clock that allocates a fixed screen-time budget to
each physical phase (see module docstring).
"""
function run_nova(params::NovaParams; fps = 30, accretion_screen_s = 8.0,
                   runaway_screen_s = 6.0, eruption_screen_s = 10.0)
    geometry = build_geometry(params)
    sol, tp = run_tnr(params)
    t_ign = find_ignition_time(sol)
    t_end = sol.t[end]

    frames = MovieFrame[]

    n_acc = max(round(Int, accretion_screen_s * fps), 2)
    acc_times = exp.(range(log(max(sol.t[1], 1.0)), log(max(t_ign, sol.t[1] + 1.0)); length = n_acc))
    for (i, t) in enumerate(acc_times)
        M_env, T, X_H, X_CNO = sol(t)
        rho, _mue, _T_F, _dr, _XHe = envelope_state([M_env, T, X_H, X_CNO], tp)
        mdot = tp.mdot_profile(t)
        L = disk_luminosity(geometry.M_wd, geometry.R_wd, mdot)
        t_movie = (i - 1) / fps
        push!(frames, MovieFrame(:accretion, t, t_movie, M_env, T, rho, X_H, X_CNO, L, mdot, geometry.R_wd))
    end

    n_run = max(round(Int, runaway_screen_s * fps), 2)
    run_times = exp.(range(log(max(t_ign, 1.0)), log(max(t_end, t_ign + 1.0)); length = n_run))
    for (i, t) in enumerate(run_times)
        M_env, T, X_H, X_CNO = sol(t)
        rho, _mue, _T_F, _dr, _XHe = envelope_state([M_env, T, X_H, X_CNO], tp)
        mdot = tp.mdot_profile(t)
        L = disk_luminosity(geometry.M_wd, geometry.R_wd, mdot)
        t_movie = accretion_screen_s + (i - 1) / fps
        push!(frames, MovieFrame(:runaway, t, t_movie, M_env, T, rho, X_H, X_CNO, L, mdot, geometry.R_wd))
    end

    M_env_peak, T_peak, X_H_peak, X_CNO_peak = sol(t_end)
    rho_peak, _mue, _T_F, _dr, _XHe = envelope_state([M_env_peak, T_peak, X_H_peak, X_CNO_peak], tp)
    eruption = run_eruption(geometry.M_wd, geometry.R_wd, T_peak, M_env_peak, rho_peak, X_H_peak;
                             n = max(round(Int, eruption_screen_s * fps), 2))
    L_peak = maximum(es.L for es in eruption)
    v_ej = eruption[1].v_ej

    for (i, es) in enumerate(eruption)
        T_eff = (es.L / (4 * pi * es.R_photosphere^2 * sigma_SB))^(1 / 4)
        t_movie = accretion_screen_s + runaway_screen_s + (i - 1) / fps
        push!(frames, MovieFrame(:eruption, t_end + es.t, t_movie, es.M_env, T_eff, rho_peak,
                                  X_H_peak, X_CNO_peak, es.L, 0.0, es.R_photosphere))
    end

    return NovaSolution(params, geometry, frames, t_ign, T_peak, L_peak, v_ej, fps)
end

total_movie_duration(sol::NovaSolution) = sol.frames[end].t_movie
