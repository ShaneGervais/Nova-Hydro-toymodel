# Drives `visualization.jl`'s scene through every frame and records it to a
# movie file via `Makie.record`. Backend-agnostic: activate `GLMakie` (GPU,
# interactive — the default for real use) or `CairoMakie` (CPU, headless —
# useful in CI/servers without a display) before calling.

"""
    make_movie(sol::NovaSolution, outfile = "nova.mp4")

Renders `sol` (from `run_nova`) to `outfile` at `sol.fps` frames per second.
Requires a Makie backend (`GLMakie`/`CairoMakie`) to already be loaded and
activated.
"""
function make_movie(sol::NovaSolution, outfile = "nova.mp4")
    fig, update! = build_scene(sol)
    n = length(sol.frames)
    Makie.record(fig, outfile, 1:n; framerate = sol.fps) do i
        update!(i)
    end
    return outfile
end
