# Physical constants and unit scales, all in CGS (cm, g, s, K, erg).
# Astrophysics codes conventionally work in CGS; keeping the numerical core in
# plain Float64 CGS avoids friction with the ODE integrator, while every public
# NovaParams field is named/commented with its physical meaning so unit errors
# are easy to catch by inspection. Unitful.jl is used in test/runtests.jl to
# independently cross-check these constants dimensionally.

const G      = 6.674_30e-8        # gravitational constant, cm^3 g^-1 s^-2
const c_light = 2.997_924_58e10    # speed of light, cm/s
const sigma_SB = 5.670_374e-5      # Stefan-Boltzmann, erg cm^-2 s^-1 K^-4
const a_rad   = 7.565_723e-15      # radiation constant, erg cm^-3 K^-4  (= 4*sigma_SB/c)
const k_B    = 1.380_649e-16       # Boltzmann constant, erg/K
const m_u    = 1.660_539e-24       # atomic mass unit, g
const m_e    = 9.109_384e-28       # electron mass, g
const hbar   = 1.054_571e-27       # reduced Planck constant, erg s
const N_A    = 6.022_140_76e23     # Avogadro's number, 1/mol

# Astronomical scales
const Msun = 1.988_47e33   # g
const Rsun = 6.957e10      # cm
const Lsun = 3.828e33      # erg/s
const yr_s  = 3.155_76e7    # Julian year, s (365.25 d)
const day_s = 8.64e4        # s

# Derived / frequently reused combinations
const G_Msun = G * Msun            # cm^3 s^-2, handy for GM_wd = m_wd_in_Msun * G_Msun
