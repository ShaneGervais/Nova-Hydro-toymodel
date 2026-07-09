module NovaSim

import OrdinaryDiffEq
import OrdinaryDiffEqRosenbrock
import SciMLBase
import Roots
import Makie
using Printf: @sprintf

include("constants.jl")
include("params.jl")
include("orbit.jl")
include("accretion.jl")
include("envelope.jl")
include("tnr.jl")
include("eruption.jl")
include("simulation.jl")
include("visualization.jl")
include("recording.jl")

export NovaParams, CompanionType, MainSequence, Subgiant, Giant,
       WDComposition, CO, ONe,
       ConstantMdot, RampMdot, FlickeringMdot,
       run_nova, NovaSolution, make_movie, build_scene,
       # physics building blocks useful for interactive exploration/teaching
       orbital_period, orbital_separation, roche_lobe_radius, l1_point,
       roche_lobe_contour, wd_radius, companion_radius, separation_for_roche_filling,
       circularization_radius, disk_forms, disk_outer_radius, disk_effective_temperature,
       stream_trajectory, run_tnr, run_eruption

end # module NovaSim
