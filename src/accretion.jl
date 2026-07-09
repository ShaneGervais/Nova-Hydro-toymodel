# Mass transfer stream, disk formation, and steady-state accretion disk
# structure.

"""
    circularization_radius(M1, M2, a) -> R_circ [cm]

Radius at which gas leaving L1 would settle into a circular Keplerian orbit
around the WD (mass M1), found from conservation of specific angular momentum:
gas at L1 corotates with the binary at distance `b1` from the WD, carrying
j = Ω b1²; equating to Keplerian j = √(G M1 R_circ) gives
`R_circ = Ω² b1⁴ / (G M1)`. See Frank, King & Raine, *Accretion Power in
Astrophysics*, §4.4. Here `b1` is computed exactly from `l1_point`/`wd_position`
rather than the usual polynomial fit, so it stays consistent with the Roche
geometry used elsewhere in this package.
"""
function circularization_radius(M1, M2, a)
    xw, _ = wd_position(M1, M2, a)
    x1 = l1_point(M1, M2, a)
    b1 = abs(x1 - xw)
    omega2 = G * (M1 + M2) / a^3
    return omega2 * b1^4 / (G * M1)
end

"""
    disk_forms(M1, M2, a, wd_composition) -> Bool

True if the circularization radius exceeds the WD radius, i.e. the stream
misses the star and must form a disk rather than impacting directly.
"""
function disk_forms(M1, M2, a, wd_composition::WDComposition)
    return circularization_radius(M1, M2, a) > wd_radius(M1, wd_composition)
end

"""
    disk_outer_radius(M1, M2, a) -> R_out [cm]

Tidal truncation caps the disk near ~0.9 of the WD's Roche lobe; viscous
spreading typically pushes the outer edge to ~1.7× the circularization radius.
We take the smaller of the two, as is standard for CV disks (Frank, King &
Raine; Warner 1995, *Cataclysmic Variable Stars*).
"""
function disk_outer_radius(M1, M2, a)
    r_circ = circularization_radius(M1, M2, a)
    r_tidal = 0.9 * roche_lobe_radius(M1, M2, a)
    return min(1.7 * r_circ, r_tidal)
end

"""
    stream_trajectory(M1, M2, a; kick=1.0e3, tmax_factor=6.0, n=400) -> (xs, ys, ts)

Ballistic path of a Lagrangian test particle released essentially at rest from
L1, integrated in the corotating binary frame under gravity from both stars
plus the centrifugal and Coriolis pseudo-forces (the restricted circular
three-body problem — the standard treatment of the mass-transfer stream, e.g.
Lubow & Shu 1975). `kick` is a tiny inward velocity perturbation [cm/s] needed
because L1 is a saddle point of the effective potential (zero velocity is an
unstable equilibrium there in principle). Integration stops when the particle
reaches the WD's Roche lobe boundary near the WD or after `tmax_factor`
orbital periods.
"""
function stream_trajectory(M1, M2, a; kick = 1.0e3, tmax_factor = 6.0, n = 400)
    xw, yw = wd_position(M1, M2, a)
    xc, yc = companion_position(M1, M2, a)
    omega2 = G * (M1 + M2) / a^3
    omega = sqrt(omega2)
    x1 = l1_point(M1, M2, a)

    function rhs!(du, u, _p, _t)
        x, y, vx, vy = u
        r1_3 = ((x - xw)^2 + (y - yw)^2)^1.5
        r2_3 = ((x - xc)^2 + (y - yc)^2)^1.5
        ax = -G * M1 * (x - xw) / r1_3 - G * M2 * (x - xc) / r2_3 + 2 * omega * vy + omega2 * x
        ay = -G * M1 * (y - yw) / r1_3 - G * M2 * (y - yc) / r2_3 - 2 * omega * vx + omega2 * y
        du[1] = vx
        du[2] = vy
        du[3] = ax
        du[4] = ay
    end

    # Initial velocity points from L1 toward the WD, i.e. down the potential
    # saddle, with a small magnitude compared to orbital velocities.
    v0x = -kick * sign(x1 - xw)
    u0 = [x1, 0.0, v0x, 0.0]
    period = 2pi / omega
    tspan = (0.0, tmax_factor * period)

    r_wd_lobe = 0.2 * roche_lobe_radius(M1, M2, a)  # "close enough to the WD" stop radius
    condition(u, _t, _integrator) = ((u[1] - xw)^2 + (u[2] - yw)^2) - r_wd_lobe^2
    affect!(integrator) = SciMLBase.terminate!(integrator)
    cb = OrdinaryDiffEq.ContinuousCallback(condition, affect!)

    prob = OrdinaryDiffEq.ODEProblem(rhs!, u0, tspan)
    sol = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.Tsit5(); callback = cb, dtmax = period / 200)

    ts = range(sol.t[1], sol.t[end]; length = n)
    xs = [sol(t)[1] for t in ts]
    ys = [sol(t)[2] for t in ts]
    return xs, ys, collect(ts)
end

"""
    disk_effective_temperature(r, M1, Mdot, R_wd) -> Teff [K]

Steady-state Shakura & Sunyaev (1973) α-disk effective temperature profile for
a disk accreting at rate `Mdot` onto a star of mass `M1` and radius `R_wd`:

    σ Teff(r)^4 = (3 G M1 Mdot)/(8π r³) · [1 - √(R_wd/r)]

Valid for R_wd ≤ r ≤ R_out. Steady-state is used rather than solving the
viscous diffusion PDE because the viscous timescale is much shorter than the
timescale over which `Ṁ(t)` varies in this model — the disk relaxes ~instantly
to the local Ṁ (documented simplification).
"""
function disk_effective_temperature(r, M1, Mdot, R_wd)
    r < R_wd && return 0.0
    flux = 3 * G * M1 * Mdot / (8 * pi * r^3) * (1 - sqrt(R_wd / r))
    return (max(flux, 0.0) / sigma_SB)^(1 / 4)
end

"""
    disk_luminosity(M1, R_wd, Mdot) -> L [erg/s]

Total accretion luminosity released in the disk down to the WD surface,
L = G M1 Ṁ / (2 R_wd) (half the total accretion energy; the other half is
released in the boundary layer where the disk meets the WD surface).
"""
disk_luminosity(M1, R_wd, Mdot) = G * M1 * Mdot / (2 * R_wd)
