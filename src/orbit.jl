# Binary orbital mechanics and Roche geometry.
#
# Convention: star 1 is always the white dwarf (accretor, mass M1), star 2 is
# the companion (donor, mass M2). Positions are in the corotating frame with
# the WD at negative x, companion at positive x, separated by `a`, center of
# mass at the origin, orbit in the x-y plane.

"""
    orbital_period(M1, M2, a) -> P [s]

Kepler's third law: P = 2π√(a³/(G(M1+M2))).
"""
orbital_period(M1, M2, a) = 2pi * sqrt(a^3 / (G * (M1 + M2)))

"""
    orbital_separation(M1, M2, P) -> a [cm]

Inverse of `orbital_period`.
"""
orbital_separation(M1, M2, P) = cbrt(G * (M1 + M2) * (P / (2pi))^2)

"""
    roche_lobe_radius(M_this, M_other, a) -> R_L [cm]

Eggleton (1983) approximation for the volume-equivalent Roche lobe radius of
the star with mass `M_this`, orbiting `M_other` at separation `a`. Accurate to
~1% for any mass ratio. Apply with `M_this` = the star whose lobe you want.
"""
function roche_lobe_radius(M_this, M_other, a)
    q = M_this / M_other
    q23 = q^(2 / 3)
    return a * 0.49 * q23 / (0.6 * q23 + log(1 + q^(1 / 3)))
end

"""
    roche_potential(x, y, M1, M2, a; frame_origin=:com)

Dimensionless-free Roche potential (per unit mass) in the corotating frame,
in erg/g, with the WD (mass M1) at (-a*M2/(M1+M2), 0) and the companion
(mass M2) at (a*M1/(M1+M2), 0), i.e. positions measured from the center of
mass, consistent with `wd_position`/`companion_position` below.

    Φ(x,y) = -G M1/r1 - G M2/r2 - ½ ω² (x² + y²)

with ω² = G(M1+M2)/a³ (Kepler), r1,r2 the distances to each star.
"""
function roche_potential(x, y, M1, M2, a)
    xw, yw = wd_position(M1, M2, a)
    xc, yc = companion_position(M1, M2, a)
    r1 = sqrt((x - xw)^2 + (y - yw)^2)
    r2 = sqrt((x - xc)^2 + (y - yc)^2)
    omega2 = G * (M1 + M2) / a^3
    return -G * M1 / r1 - G * M2 / r2 - 0.5 * omega2 * (x^2 + y^2)
end

wd_position(M1, M2, a) = (-a * M2 / (M1 + M2), 0.0)
companion_position(M1, M2, a) = (a * M1 / (M1 + M2), 0.0)

"""
    l1_point(M1, M2, a) -> x [cm]

x-coordinate (on the line joining the stars) of the inner Lagrange point L1,
found by root-finding the effective 1D force balance along the line of
centers, using `Roots.jl`. Returns a coordinate in the same center-of-mass
frame as `wd_position`/`companion_position`.
"""
function l1_point(M1, M2, a)
    xw, _ = wd_position(M1, M2, a)
    xc, _ = companion_position(M1, M2, a)
    omega2 = G * (M1 + M2) / a^3
    # Explicit derivative of Φ(x,0) w.r.t. x (careful with signs of r):
    #  Φ(x,0) = -GM1/|x-xw| - GM2/|x-xc| - ½ ω² x²
    #  dΦ/dx  =  GM1*sign(x-xw)/(x-xw)^2 + GM2*sign(x-xc)/(x-xc)^2 - ω² x
    force(x) = G * M1 * sign(x - xw) / (x - xw)^2 +
               G * M2 * sign(x - xc) / (x - xc)^2 - omega2 * x
    lo = xw + 1.0e-6 * a
    hi = xc - 1.0e-6 * a
    return Roots.find_zero(force, (lo, hi), Roots.Brent())
end

"""
    roche_lobe_contour(M1, M2, a; n=200) -> (xs, ys)

Points on the Roche-lobe equipotential surface (the surface through L1) in the
orbital plane, traced by root-finding radial distance from each star's center
at a set of angles, for plotting the classic "figure-eight" lobe boundary.
Returns two lobes concatenated with a `NaN` separator (ready to pass straight
to `Makie.lines!`).
"""
function roche_lobe_contour(M1, M2, a; n = 200)
    x1 = l1_point(M1, M2, a)
    phi_L1 = roche_potential(x1, 0.0, M1, M2, a)
    xw, yw = wd_position(M1, M2, a)
    xc, yc = companion_position(M1, M2, a)

    function lobe_points(cx, cy, other_cx)
        xs = Float64[]
        ys = Float64[]
        r_guess = abs(x1 - cx)
        for theta in range(0, 2pi; length = n)
            dirx, diry = cos(theta), sin(theta)
            f(r) = roche_potential(cx + r * dirx, cy + r * diry, M1, M2, a) - phi_L1
            # bracket outward from a small radius until f changes sign or we
            # exceed a generous bound (guards against the open (companion-far)
            # side of a very asymmetric lobe).
            r_lo, r_hi = 1.0e-8 * a, 5.0 * abs(other_cx - cx)
            r = try
                Roots.find_zero(f, (r_lo, r_hi), Roots.Brent())
            catch
                r_guess
            end
            push!(xs, cx + r * dirx)
            push!(ys, cy + r * diry)
        end
        return xs, ys
    end

    xs1, ys1 = lobe_points(xw, yw, xc)
    xs2, ys2 = lobe_points(xc, yc, xw)
    xs = vcat(xs1, NaN, xs2)
    ys = vcat(ys1, NaN, ys2)
    return xs, ys
end

# --- Mass-radius relations for the companion (donor) -------------------------

"""
    companion_radius(M2, companion_type) -> R [cm]

Zero-age main-sequence / evolved mass-radius relations used only to *check*
Roche-lobe filling and to size the companion for visualization; not evolved in
time. `MainSequence` uses the standard low-mass power law R ∝ M^0.8 (valid
below ~1 M_sun, appropriate for CV donors). `Subgiant`/`Giant` use generously
inflated radii typical of the symbiotic/recurrent-nova donor regime.
"""
function companion_radius(M2, companion_type::CompanionType)
    m = M2 / Msun
    if companion_type == MainSequence
        return (m < 1.0 ? m^0.8 : m^0.57) * Rsun
    elseif companion_type == Subgiant
        return 3.0 * m^0.8 * Rsun
    else # Giant
        return 25.0 * m^0.6 * Rsun
    end
end

"""
    wd_radius(M1, composition) -> R [cm]

Nauenberg (1972) zero-temperature degenerate mass-radius relation:
    R = R0 * [(Mch/M)^(2/3) - (M/Mch)^(2/3)]^(1/2)
with composition-dependent (R0, Mch). CO and ONe WDs share the same
qualitative form; ONe uses a slightly higher mean molecular weight per
electron, giving a smaller radius at fixed mass.
"""
function wd_radius(M1, composition::WDComposition)
    Mch = (composition == CO ? 1.44 : 1.38) * Msun
    R0 = (composition == CO ? 0.0126 : 0.0115) * Rsun
    m = min(M1 / Mch, 0.999)
    return R0 * sqrt(m^(-2 / 3) - m^(2 / 3))
end

"""
    separation_for_roche_filling(M1, M2, companion_type) -> a [cm]

Solves for the orbital separation at which the companion (mass M2, type
`companion_type`) exactly fills its Roche lobe, i.e.
`companion_radius(M2, companion_type) == roche_lobe_radius(M2, M1, a)`. This is
the natural default separation for a mass-transferring binary.
"""
function separation_for_roche_filling(M1, M2, companion_type::CompanionType)
    R2 = companion_radius(M2, companion_type)
    q = M2 / M1
    q23 = q^(2 / 3)
    f = 0.49 * q23 / (0.6 * q23 + log(1 + q^(1 / 3)))
    return R2 / f
end
