using FunWithMIPS
#=using ArgParse

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--phi"
            arg_type = Float64
            default = 0.5
        "--"
            arg_type = Int64
            default = 0
        "--traj_filename"
            arg_type = String
            default = "TEST_OUTPUTS/trajectories_tetramer_test.txt"
    end

    return parse_args(s)
end=#

function run!()
    radius = 0.5

    phi = 0.6
    L = 50
    f_propulsion = 0.5
    D_rot = 0.001
    kT = 0.0

    num_steps = 100000
    dt = 0.01

    traj_partial_filename = "TEST_OUTPUTS/trajectory_partial.out"
    animation_partial_filename = "TEST_OUTPUTS/animation_partial.mp4"
    traj_full_filename = "TEST_OUTPUTS/trajectory_full.out"
    animation_full_filename = "TEST_OUTPUTS/animation_full.mp4"
    num_frames_to_save = 300
    save_interval_partial = trunc(Int64, 0.5 / (f_propulsion * dt))
    save_start_partial = num_steps - save_interval_partial * num_frames_to_save
    save_interval_full = trunc(Int64, num_steps / num_frames_to_save)
    save_start_full = num_steps - save_interval_full * num_frames_to_save

    lattice_spacing = sqrt(2 * pi * radius^2 / (sqrt(3) * phi))
    Nx = trunc(Int64, L / lattice_spacing)
    initial_positions, box = initialize_triangular_lattice(lattice_spacing, (Nx, trunc(Int64, Nx / sqrt(3))))
    particles = Vector{Particle}()
    for position in initial_positions
        push!(particles, Particle(position, ActiveBrownian(f_propulsion, D_rot)))
    end
    println("Number of particles: $(length(particles))")
    println("Packing fraction: $(phi)")
    println("Box: $(box)")

    cell_list, padding = CellList(particles, 2 * radius, box; approx_padding = 0.25)
    cell_list.update_interval = maximum([1, floor(Int64, 0.90 * padding / (f_propulsion * dt))])
    println("CellList update interval set to $(cell_list.update_interval)")
    k = 10.0
    lj = HarmonicRepulsion(k, 2 * radius, cell_list, box)

    brownian = BrownianDynamics(particles, dt, kT, box)

    system = System(particles, [lj], [cell_list], brownian)
    partial_trajectory = TrajectoryContainer(save_interval_partial, box; start = save_start_partial)
    full_trajectory = TrajectoryContainer(save_interval_full, box; start = save_start_full)
    run_simulation!(system, [partial_trajectory, full_trajectory], num_steps)
    save!(partial_trajectory, traj_partial_filename)
    save!(full_trajectory, traj_full_filename)

    #visualize!(partial_trajectory)
    #visualize!(full_trajectory)
    save_movie!(partial_trajectory, 30, animation_partial_filename)
    save_movie!(full_trajectory, 30, animation_full_filename)
end
run!()