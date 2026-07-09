# Multi-panel movie scene, written against the backend-agnostic `Makie` API.
# The caller picks a concrete backend (`using GLMakie` for interactive/GPU
# rendering, `using CairoMakie` for headless CPU rendering) and activates it
# before calling `build_scene`/`make_movie` — this file never imports a
# specific backend, so both work unchanged (Makie's intended design).
#
# Layout:
#   [ orbital view (binary, Roche lobes, stream, disk/ejecta) | envelope cutaway ]
#   [        luminosity(t)        |   temperature(t)   |    Ṁ(t)              ]
# with a clock readout (phase name, physical time elapsed, movie time) drawn
# over the orbital panel.

"""
    format_duration(seconds) -> String

Human-readable physical-time label spanning seconds to millennia, used for the
on-screen clock readout.
"""
function format_duration(seconds)
    s = abs(seconds)
    if s < 1
        return @sprintf("%.1f ms", seconds * 1000)
    elseif s < 60
        return @sprintf("%.1f s", seconds)
    elseif s < 3600
        return @sprintf("%.1f min", seconds / 60)
    elseif s < day_s
        return @sprintf("%.1f hr", seconds / 3600)
    elseif s < yr_s
        return @sprintf("%.1f d", seconds / day_s)
    else
        return @sprintf("%.2f yr", seconds / yr_s)
    end
end

phase_label(phase::Symbol) = phase == :accretion ? "Accretion (Roche-lobe overflow)" :
                              phase == :runaway ? "Thermonuclear runaway" : "Eruption / ejection"

"""
    zone_temperatures(T_base, n_zones) -> Vector{Float64}

Illustrative radial temperature profile for the envelope cutaway: an
exponential decline from the hot burning base outward. Envelope thickness is
exaggerated ~1000x for visibility — real nova envelopes are a geometrically
thin skin (~10^-4 R_wd) that would be invisible drawn to scale (documented
artistic license, see README).
"""
zone_temperatures(T_base, n_zones) = [T_base * exp(-3 * (i - 1) / max(n_zones - 1, 1)) for i in 1:n_zones]

const DISK_N_POINTS = 2000

"""
    build_scene(sol::NovaSolution) -> (fig, update!)

Builds the figure and all plot elements for `sol`, returning the `Makie`
figure and an `update!(i::Int)` closure that mutates every Observable to
frame `i`. `recording.jl` drives `update!` inside `Makie.record`.
"""
function build_scene(sol::NovaSolution)
    geom = sol.geometry
    frames = sol.frames

    fig = Makie.Figure(size = (1600, 950), backgroundcolor = :black)

    # --- Orbital panel -----------------------------------------------------
    # `span` is the wide, "establishing shot" framing (both stars + Roche
    # lobes + stream). The camera dynamically zooms in toward the WD during
    # the runaway (where all the interesting physics is, but which is tiny
    # next to the orbital scale) and pulls back out during the eruption to
    # follow the expanding ejecta — see `camera_span` below. Axis limits are
    # therefore *not* fixed at construction; `update!` drives them every frame.
    span = 1.3 * max(abs(geom.wd_pos[1]), abs(geom.companion_pos[1]), geom.disk_r_out, geom.R_companion)
    ax_orbit = Makie.Axis(fig[1, 1]; backgroundcolor = :black, aspect = Makie.DataAspect(),
                           title = "Binary system", titlecolor = :white,
                           xticksvisible = false, yticksvisible = false,
                           xticklabelsvisible = false, yticklabelsvisible = false)

    Makie.lines!(ax_orbit, geom.roche_xs, geom.roche_ys; color = (:dodgerblue, 0.6), linewidth = 1.5)
    Makie.lines!(ax_orbit, geom.stream_xs, geom.stream_ys; color = (:orange, 0.7), linewidth = 1.5)

    theta_disk = 2pi .* rand(DISK_N_POINTS)
    r_disk = geom.disk_r_in .+ (geom.disk_r_out - geom.disk_r_in) .* rand(DISK_N_POINTS) .^ 0.5
    # The WD and its disk are physically tiny next to the orbital separation
    # (R_wd/a is typically <1%) — plotted to true scale the disk would be
    # invisible. Positions are exaggerated so the disk always spans ~18% of
    # the panel; `disk_effective_temperature` below still uses the real
    # physical `r_disk`, so the temperature/color mapping stays accurate —
    # only the on-screen position is stretched.
    disk_visual_scale = 0.18 * span / geom.disk_r_out
    r_disk_visual = r_disk .* disk_visual_scale
    disk_x = Makie.Observable(geom.wd_pos[1] .+ r_disk_visual .* cos.(theta_disk))
    disk_y = Makie.Observable(geom.wd_pos[2] .+ r_disk_visual .* sin.(theta_disk))
    disk_color = Makie.Observable(fill(1.0e4, DISK_N_POINTS))
    if geom.disk_forms
        # Log color scale: disk Teff spans a couple of orders of magnitude
        # across the Ṁ range this package supports (few×10^3 K at low Ṁ up
        # to ~10^5-10^6 K for strong accretors), which a linear colorrange
        # would crush to near-black at the low end.
        Makie.scatter!(ax_orbit, disk_x, disk_y; color = disk_color, colormap = :inferno,
                        colorscale = log10, colorrange = (2.0e3, 2.0e5), markersize = 5, alpha = 0.9)
    end

    # Ejecta is drawn as an expanding, fading *ring* (not a filled disk): at
    # ~8000 km/s over the eruption's multi-day duration, the physical
    # photosphere radius quickly exceeds the whole orbital view — a filled
    # disk would just flood the panel solid. The displayed radius is capped
    # at the visible span and the ring fades out as it "exits the frame",
    # which reads as an expanding shockwave rather than a screen full of color.
    ejecta_radius = Makie.Observable(geom.R_wd)
    ejecta_alpha = Makie.Observable(0.0)
    ejecta_color = Makie.Observable(:white)
    Makie.poly!(ax_orbit, Makie.lift(r -> Makie.Circle(Makie.Point2f(geom.wd_pos...), Float32(r)), ejecta_radius);
                color = :transparent, strokecolor = Makie.lift((c, a) -> (c, a), ejecta_color, ejecta_alpha),
                strokewidth = 3)

    # WD body: drawn as a real data-space object (not a fixed-pixel marker)
    # so the dynamic camera above actually reveals it — sized off the same
    # `disk_visual_scale` exaggeration as the disk so it sits naturally at
    # the disk's inner edge. Layered glow + a surface colored by the *actual*
    # current envelope temperature (an Observable) gives the close-up shot
    # during runaway real physical content instead of a static white dot.
    wd_r = geom.R_wd * disk_visual_scale
    wd_circle(mult) = Makie.Circle(Makie.Point2f(geom.wd_pos...), Float32(wd_r * mult))
    Makie.poly!(ax_orbit, wd_circle(3.5); color = (:aliceblue, 0.10), strokewidth = 0)
    Makie.poly!(ax_orbit, wd_circle(2.0); color = (:lightskyblue, 0.20), strokewidth = 0)
    wd_surface_color = Makie.Observable(1.0e7)
    Makie.poly!(ax_orbit, wd_circle(1.0); color = wd_surface_color, colormap = :inferno,
                colorrange = (1.0e7, 4.0e8), strokewidth = 1, strokecolor = :lightblue)
    Makie.poly!(ax_orbit, wd_circle(0.4); color = :white, strokewidth = 0)

    # Companion: rendered at true physical scale, so it visibly fills its own
    # Roche lobe (touching the lobe boundary already drawn above) — a
    # physically meaningful detail, not just decoration. Three nested
    # circles approximate limb darkening (brighter toward the center).
    comp_circle(mult) = Makie.Circle(Makie.Point2f(geom.companion_pos...), Float32(geom.R_companion * mult))
    Makie.poly!(ax_orbit, comp_circle(1.0); color = :orangered, strokewidth = 0)
    Makie.poly!(ax_orbit, comp_circle(0.72); color = :orange, strokewidth = 0)
    Makie.poly!(ax_orbit, comp_circle(0.38); color = :navajowhite, strokewidth = 0)

    clock_text = Makie.Observable("")
    Makie.text!(ax_orbit, 0.02, 0.98; text = clock_text, space = :relative, color = :white,
                fontsize = 16, align = (:left, :top))

    # --- Envelope cutaway panel ---------------------------------------------
    ax_env = Makie.Axis(fig[1, 2]; backgroundcolor = :black, aspect = Makie.DataAspect(),
                         title = "WD envelope (schematic, not to radial scale)", titlecolor = :white,
                         xticksvisible = false, yticksvisible = false,
                         xticklabelsvisible = false, yticklabelsvisible = false)
    n_zones = sol.params.n_zones
    env_span = geom.R_wd * 1.5
    Makie.limits!(ax_env, -env_span, env_span, -env_span, env_span)
    zone_radii = [Makie.Observable(geom.R_wd) for _ in 1:n_zones]
    zone_colors = [Makie.Observable(1.0e7) for _ in 1:n_zones]
    for i in n_zones:-1:1
        Makie.poly!(ax_env, Makie.lift(r -> Makie.Circle(Makie.Point2f(0, 0), Float32(r)), zone_radii[i]);
                    color = zone_colors[i], colormap = :inferno, colorrange = (1.0e7, 4.0e8), strokewidth = 0)
    end
    Makie.poly!(ax_env, Makie.Circle(Makie.Point2f(0, 0), Float32(geom.R_wd * 0.3)); color = :gray10)

    # --- Strip charts --------------------------------------------------------
    t_movie_all = [f.t_movie for f in frames]
    L_all = max.([f.L for f in frames], 1.0)
    T_all = max.([f.T for f in frames], 1.0)
    mdot_all = [f.mdot * yr_s / Msun for f in frames]  # Msun/yr, easier to read

    ax_L = Makie.Axis(fig[2, 1]; backgroundcolor = :black, yscale = log10, title = "Luminosity [erg/s]",
                       titlecolor = :white, xlabel = "movie time [s]", xlabelcolor = :white,
                       ylabelcolor = :white, xticklabelcolor = :white, yticklabelcolor = :white)
    Makie.lines!(ax_L, t_movie_all, L_all; color = :yellow)
    vline_L = Makie.vlines!(ax_L, [0.0]; color = :white, linestyle = :dash)

    ax_T = Makie.Axis(fig[2, 2]; backgroundcolor = :black, yscale = log10, title = "Envelope T [K]",
                       titlecolor = :white, xlabel = "movie time [s]", xlabelcolor = :white,
                       ylabelcolor = :white, xticklabelcolor = :white, yticklabelcolor = :white)
    Makie.lines!(ax_T, t_movie_all, T_all; color = :orangered)
    vline_T = Makie.vlines!(ax_T, [0.0]; color = :white, linestyle = :dash)

    ax_M = Makie.Axis(fig[2, 3]; backgroundcolor = :black, title = "Ṁ [M_sun/yr]",
                       titlecolor = :white, xlabel = "movie time [s]", xlabelcolor = :white,
                       ylabelcolor = :white, xticklabelcolor = :white, yticklabelcolor = :white)
    Makie.lines!(ax_M, t_movie_all, mdot_all; color = :dodgerblue)
    vline_M = Makie.vlines!(ax_M, [0.0]; color = :white, linestyle = :dash)

    Makie.Label(fig[0, :], "NovaSim — classical nova cycle"; color = :white, fontsize = 22)

    # --- Dynamic camera ------------------------------------------------------
    # `span` (the wide establishing shot) holds for accretion; the camera
    # eases in toward the WD over the first ~30% of the runaway phase and
    # holds there, then eases out to a wide pull-back over the first ~40% of
    # the eruption phase — timed to arrive by the time the ejecta would
    # otherwise have already expanded past the frame, so the camera reveals
    # the expansion instead of the ring appearing to vanish.
    span_close = max(4.0 * geom.disk_r_out, 8.0 * geom.R_wd)
    span_far = 2.5 * span
    run_ts = [f.t_movie for f in frames if f.phase == :runaway]
    erupt_ts = [f.t_movie for f in frames if f.phase == :eruption]
    t_run_start, t_run_end = run_ts[1], run_ts[end]
    t_erupt_start, t_erupt_end = erupt_ts[1], erupt_ts[end]
    ease(x) = x <= 0 ? 0.0 : x >= 1 ? 1.0 : x^2 * (3 - 2x)

    function camera_span(f)
        if f.phase == :accretion
            return span
        elseif f.phase == :runaway
            prog = (f.t_movie - t_run_start) / max(t_run_end - t_run_start, 1.0e-9)
            return span + (span_close - span) * ease(clamp(prog / 0.3, 0.0, 1.0))
        else
            prog = (f.t_movie - t_erupt_start) / max(t_erupt_end - t_erupt_start, 1.0e-9)
            return span_close + (span_far - span_close) * ease(clamp(prog / 0.4, 0.0, 1.0))
        end
    end

    function update!(i::Int)
        f = frames[i]
        current_span = camera_span(f)
        Makie.xlims!(ax_orbit, -current_span, current_span)
        Makie.ylims!(ax_orbit, -current_span, current_span)

        if f.phase == :accretion || f.phase == :runaway
            mdot_now = f.mdot
            disk_color[] = disk_effective_temperature.(r_disk, geom.M_wd, mdot_now, geom.R_wd)
            ejecta_alpha[] = 0.0
        else
            disk_color[] = fill(1.0e3, DISK_N_POINTS)
            fade = f.R_photosphere <= current_span ? 1.0 : max(0.0, 1.0 - (f.R_photosphere / current_span - 1.0))
            ejecta_radius[] = min(f.R_photosphere, 1.3 * current_span)
            ejecta_alpha[] = 0.9 * fade
            ejecta_color[] = f.T > 1.0e4 ? :white : :orangered
        end
        wd_surface_color[] = f.T

        zt = zone_temperatures(f.T, n_zones)
        for i2 in 1:n_zones
            zone_radii[i2][] = geom.R_wd * (1 + 0.5 * i2 / n_zones)
            zone_colors[i2][] = zt[i2]
        end

        clock_text[] = phase_label(f.phase) * "\n" *
                       "elapsed: " * format_duration(f.t_phys) * "\n" *
                       @sprintf("T_base = %.2e K   L = %.2e erg/s", f.T, max(f.L, 1.0))

        vline_L[1] = [f.t_movie]
        vline_T[1] = [f.t_movie]
        vline_M[1] = [f.t_movie]
        return nothing
    end

    update!(1)
    return fig, update!
end
