# The "common" classical nova default: a 0.8 M_sun CO white dwarf accreting
# at a modest steady rate from a 0.4 M_sun main-sequence donor filling its
# Roche lobe. Produces a single movie of the full accretion -> thermonuclear
# runaway -> eruption cycle.
#
# GLMakie needs a real display/GPU (it will not even load headlessly, since
# its GLFW dependency initializes a window context on `using`). On a machine
# with a display this is the right backend for interactive/GPU-accelerated
# rendering. For headless servers/CI, swap the next line for `using CairoMakie`
# — every other line of this script, and of NovaSim itself, is unchanged,
# since the plotting code in src/visualization.jl is written against the
# backend-agnostic `Makie` API.
using GLMakie

using NovaSim

params = NovaParams()  # all fields have physically-motivated defaults; try e.g.
                        # NovaParams(M_wd = 1.0Msun, mixing_fraction = 0.5) for
                        # a more massive, more strongly-mixed ONe nova instead.

sol = run_nova(params)

@info "Ignition at" t_ignition_yr = sol.t_ignition / NovaSim.yr_s
@info "Peak temperature [K]" sol.T_peak
@info "Peak luminosity [erg/s]" sol.L_peak
@info "Ejecta velocity [km/s]" v_ejecta_kms = sol.v_ejecta / 1.0e5

make_movie(sol, joinpath(@__DIR__, "classical_nova_co_wd.mp4"))
