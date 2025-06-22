using CSV
using DataFrames
using Plots
using Printf

function circleShape(xc, yc, rx, ry; N=100)
    θ = range(0, 2π, length=N)
    x = xc .+ rx .* cos.(θ)
    y = yc .+ ry .* sin.(θ)
    return x, y
end

function animate_particles(df::DataFrame, output_filename::String; plots_to_show::Vector{Symbol} = [:position, :velocity, :width], frame_step::Int = 1)
    total_cols = size(df, 2)
    n_particles = (total_cols - 2) ÷ 3

    pos_cols = 1:n_particles
    vel_cols = (n_particles+1):(2n_particles)
    wid_cols = (2n_particles+1):(3n_particles)
    time_col = total_cols - 1
    dt_col = total_cols

    initial_widths = [df[1, i] for i in wid_cols]
    initial_radii = initial_widths ./ 2
    max_r = maximum(initial_radii)

    time = df[:, time_col]
    h = df[2, time_col] - df[1, time_col]

    # Determine full axis limits for consistent scaling
    pos_min = minimum([minimum(df[:, j]) for j in pos_cols]) - 1.0
    pos_max = maximum([maximum(df[:, j]) for j in pos_cols]) + 1.0
    vel_min = minimum([minimum(df[:, j]) for j in vel_cols])
    vel_max = maximum([maximum(df[:, j]) for j in vel_cols])

    # Precompute deformation (initial_width - width)
    deformation_df = [initial_widths[j] .- df[:, wid_cols[j]] for j in 1:n_particles]
    deform_min = minimum(reduce(vcat, deformation_df))
    deform_max = maximum(reduce(vcat, deformation_df))

    anim = @animate for i in 1:frame_step:size(df, 1)
        positions = [df[i, j] for j in pos_cols]
        velocities = [df[i, j] for j in vel_cols]
        widths = [df[i, j] for j in wid_cols]
        radii_x = widths ./ 2
        radii_y = initial_radii

        # Particle outlines
        Xc, Yc = Vector{Vector{Float64}}(), Vector{Vector{Float64}}()
        for j in 1:n_particles
            x_circle, y_circle = circleShape(positions[j], max_r, radii_x[j], radii_y[j])
            push!(Xc, x_circle)
            push!(Yc, y_circle)
        end

        # Particle colors by deformation
        colors = [widths[j] < initial_widths[j] ? "red" : "blue" for j in 1:n_particles]

        Xs = [[positions[j]] for j in 1:n_particles]
        Ys = [[max_r] for _ in 1:n_particles]

        # Particle plot
        p_particles = plot(xlim = (pos_min,pos_max),
                           ylim = (0, 2*max_r),
                           aspect_ratio = :equal,
                           legend = false)
        for j in 1:n_particles
            plot!(p_particles, Xc[j], Yc[j],
                  linecolor = colors[j],
                  linewidth = 3,
                  label = false)
        end
        title!(p_particles, @sprintf("t = %.2f", df[i, time_col]))

        # Time-series plots
        p_pos = plot(title = "Positions", xlim = (minimum(time), maximum(time)), ylim = (pos_min, pos_max), legend = false)
        for j in pos_cols
            plot!(p_pos, time[1:i], df[1:i, j], label = "", lw = 3)
        end

        p_vel = plot(title = "Velocities", xlim = (minimum(time), maximum(time)), ylim = (vel_min, vel_max), legend = false)
        for j in vel_cols
            plot!(p_vel, time[1:i], df[1:i, j], label = "", lw = 3)
        end

        p_def = plot(title = "Deformation (initial - current width)",
                     xlim = (minimum(time), maximum(time)),
                     ylim = (deform_min, deform_max),
                     legend = false)
        for j in 1:n_particles
            plot!(p_def, time[1:i], deformation_df[j][1:i], label = "", lw = 3)
        end

        l = @layout [grid(1,3); b{0.5h}]
        plot(p_pos, p_vel, p_def, p_particles, layout = l, size = (1200, 800))

        # Progress bar
        progress = i / size(df, 1) * 100
        print("\rAnimating frame $(i)/$(size(df,1)) — $(round(progress, digits=1))%")
        flush(stdout)
    end

    println("\nSaving animation to $(output_filename)...")
    gif(anim, output_filename, fps = 20)
    println("Animation saved.")
end
