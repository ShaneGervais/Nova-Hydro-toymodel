# NovaParams: the single struct carrying every physical knob the simulation
# exposes. All fields are plain Float64 in CGS units (see constants.jl) unless
# noted; use the `Msun`, `Rsun`, `yr_s`, etc. scale constants at the call site
# to convert, e.g. `M_wd = 0.8Msun`.

"""
    CompanionType

Donor-star type, selects the mass-radius relation used to check Roche-lobe
filling (`radius_main_sequence`, `radius_subgiant`, `radius_giant` in orbit.jl).
"""
@enum CompanionType MainSequence Subgiant Giant

"""
    WDComposition

White dwarf core composition. Selects mass-radius normalization, allowed mass
range, and the simplified nuclear energy-generation / opacity coefficients used
during the thermonuclear runaway (tnr.jl). CO is the common classical-nova case;
ONe ("neon nova") WDs are more massive and burn hotter.
"""
@enum WDComposition CO ONe

"""
    MdotProfile

A callable object giving the donor->WD mass transfer rate Ṁ(t) [g/s] as a
function of physical time t [s] since the simulation start. Use one of the
constructors below, or pass any `t -> value` function directly.
"""
abstract type MdotProfile end

struct ConstantMdot <: MdotProfile
    mdot::Float64  # g/s
end
(p::ConstantMdot)(t) = p.mdot

"""
    RampMdot(mdot0, rate)

Linearly increasing (or decreasing, for negative `rate`) accretion rate:
Ṁ(t) = mdot0 + rate*t, floored at zero.
"""
struct RampMdot <: MdotProfile
    mdot0::Float64  # g/s at t=0
    rate::Float64   # g/s per s
end
(p::RampMdot)(t) = max(p.mdot0 + p.rate * t, 0.0)

"""
    FlickeringMdot(mean, amplitude, period; phase=0.0)

Oscillating accretion rate representing disk-instability / flickering
variability: Ṁ(t) = mean + amplitude*sin(2π t/period + phase), floored at zero.
"""
struct FlickeringMdot <: MdotProfile
    mean::Float64
    amplitude::Float64
    period::Float64  # s
    phase::Float64
end
FlickeringMdot(mean, amplitude, period; phase = 0.0) = FlickeringMdot(mean, amplitude, period, phase)
(p::FlickeringMdot)(t) = max(p.mean + p.amplitude * sin(2pi * t / p.period + p.phase), 0.0)

# Allow bare functions too.
const AnyMdotProfile = Union{MdotProfile, Function}

"""
    NovaParams(; kwargs...)

Every physical parameter of a nova simulation. Defaults describe a common
classical nova: a 0.8 M_sun CO white dwarf accreting at a modest steady rate
from a 0.4 M_sun main-sequence donor filling its Roche lobe.

Fields
- `M_wd`             : white dwarf mass [g]
- `wd_composition`    : `CO` or `ONe`
- `M_companion`       : donor mass [g]
- `companion_type`    : `MainSequence`, `Subgiant`, or `Giant`
- `mdot_profile`      : `MdotProfile` (or `t -> g/s` function) giving Ṁ(t)
- `separation`        : orbital separation [cm]; `NaN` => derived from Roche-lobe filling
- `mixing_fraction`   : fraction (0-1) of core WD material dredged into the burning
                         envelope by convection during runaway
- `convection_alpha`  : convective mixing-efficiency parameter (mixing-length-like,
                         order-unity; larger = faster homogenization)
- `base_density`      : envelope base density override [g/cm^3]; `NaN` => derived
                         from hydrostatic equilibrium at the WD surface
- `core_temperature`  : WD core/isothermal temperature before accretion [K]
- `n_zones`           : number of Lagrangian shells in the envelope TNR model
"""
@kwdef struct NovaParams
    M_wd::Float64            = 0.8 * Msun
    wd_composition::WDComposition = CO
    M_companion::Float64     = 0.4 * Msun
    companion_type::CompanionType = MainSequence
    mdot_profile::AnyMdotProfile = ConstantMdot(5.0e-10 * Msun / yr_s)
    separation::Float64       = NaN
    mixing_fraction::Float64  = 0.25
    convection_alpha::Float64 = 1.5
    base_density::Float64     = NaN
    core_temperature::Float64 = 1.0e7
    n_zones::Int               = 8
end

mdot(p::NovaParams, t) = p.mdot_profile(t)
