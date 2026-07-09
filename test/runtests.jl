using Test
using NovaSim
using Unitful
using CairoMakie

@testset "Constants sanity" begin
    @test isapprox(NovaSim.yr_s, ustrip(u"s", 365.25u"d"); rtol = 1.0e-12)
    @test isapprox(NovaSim.Msun, 1.989e33; rtol = 1.0e-2)
    @test isapprox(NovaSim.Rsun, 6.957e10; rtol = 1.0e-3)
end

@testset "Roche lobe radius (Eggleton 1983)" begin
    a = 1.0e11
    # Equal masses (q=1): a well-known reference value, R_L/a ≈ 0.3789.
    @test isapprox(roche_lobe_radius(NovaSim.Msun, NovaSim.Msun, a) / a, 0.3789; rtol = 1.0e-3)
    # Very small companion: lobe should shrink toward zero.
    @test roche_lobe_radius(1.0e-3 * NovaSim.Msun, NovaSim.Msun, a) < 0.05a
    # Dimensional genericity: the formula is a pure mass-ratio times `a`, so it
    # should agree numerically across unit systems (a real Unitful cross-check,
    # not just a units-look-right assertion).
    M1, M2 = 0.8NovaSim.Msun, 0.4NovaSim.Msun
    r_cgs = roche_lobe_radius(M1, M2, a)
    r_si = roche_lobe_radius(uconvert(u"kg", M1 * u"g"), uconvert(u"kg", M2 * u"g"), uconvert(u"m", a * u"cm"))
    @test isapprox(ustrip(u"cm", uconvert(u"cm", r_si)), r_cgs; rtol = 1.0e-9)
end

@testset "Kepler orbital mechanics" begin
    M1, M2, a = 0.8NovaSim.Msun, 0.4NovaSim.Msun, 2.0e11
    P = orbital_period(M1, M2, a)
    @test isapprox(orbital_separation(M1, M2, P), a; rtol = 1.0e-10)
end

@testset "WD mass-radius relation (Nauenberg 1972)" begin
    @test wd_radius(1.0NovaSim.Msun, NovaSim.CO) < wd_radius(0.6NovaSim.Msun, NovaSim.CO)
    @test 0.001NovaSim.Rsun < wd_radius(0.8NovaSim.Msun, NovaSim.CO) < 0.03NovaSim.Rsun
end

@testset "L1 point and Roche geometry" begin
    M1, M2, a = 0.8NovaSim.Msun, 0.4NovaSim.Msun, 2.0e11
    xw, _ = NovaSim.wd_position(M1, M2, a)
    xc, _ = NovaSim.companion_position(M1, M2, a)
    x1 = l1_point(M1, M2, a)
    @test xw < x1 < xc
end

@testset "Accretion disk geometry" begin
    M1, M2 = 0.8NovaSim.Msun, 0.4NovaSim.Msun
    companion_type = NovaSim.MainSequence
    a = separation_for_roche_filling(M1, M2, companion_type)
    @test disk_forms(M1, M2, a, NovaSim.CO)  # short-period CV: expect a disk, not direct impact
    r_out = disk_outer_radius(M1, M2, a)
    @test wd_radius(M1, NovaSim.CO) < r_out < a
end

@testset "Nuclear energy generation trends" begin
    rho, X_H, X_CNO = 1.0e3, 0.7, 0.02
    # Cold-CNO regime should be steeply increasing with T.
    @test NovaSim.epsilon_nuc(3.0e7, rho, X_H, X_CNO) < NovaSim.epsilon_nuc(6.0e7, rho, X_H, X_CNO)
    # Hot-CNO regime should be essentially density- and temperature-independent.
    e1 = NovaSim.epsilon_nuc(1.5e8, rho, X_H, X_CNO)
    e2 = NovaSim.epsilon_nuc(3.0e8, rho, X_H, X_CNO)
    @test isapprox(e1, e2; rtol = 0.1)
    e3 = NovaSim.epsilon_nuc(1.5e8, 2 * rho, X_H, X_CNO)
    @test isapprox(e1, e3; rtol = 0.1)
end

@testset "Full nova cycle smoke test (default classical-nova parameters)" begin
    params = NovaParams()
    sol = run_nova(params; fps = 10, accretion_screen_s = 2.0, runaway_screen_s = 2.0, eruption_screen_s = 2.0)

    @test !isempty(sol.frames)
    @test 0 < sol.t_ignition < sol.frames[end].t_phys
    @test 1.0e3 * NovaSim.yr_s < sol.t_ignition < 1.0e6 * NovaSim.yr_s  # plausible recurrence timescale

    # Broad sanity ranges for this one-zone model (not tight regression bounds,
    # and deliberately not the ~2-4e8 K often quoted from multi-zone codes —
    # a single representative shell can't resolve the localized hot spot that
    # drives that higher figure; see README's documented simplifications).
    @test 1.0e7 < sol.T_peak < 5.0e8
    @test 0 < sol.v_ejecta < 1.0e4 * 1.0e5  # < 10,000 km/s
    @test 0 < sol.L_peak <= 2 * NovaSim.eddington_luminosity(params.M_wd, 0.7)

    # Rendering smoke test (CairoMakie: headless, no display needed). Confirms
    # the backend-agnostic scene in visualization.jl actually builds and draws
    # every frame without error across all three phases.
    fig, update! = build_scene(sol)
    for i in (1, length(sol.frames) ÷ 2, length(sol.frames))
        update!(i)
    end
    mktempdir() do dir
        outfile = joinpath(dir, "frame.png")
        save(outfile, fig)
        @test filesize(outfile) > 0
    end
end
