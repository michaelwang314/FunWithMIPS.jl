using FunWithMIPS
#=using ArgParse

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--phi"
            arg_type = Float64
            default = 0.5
        "--num_steps"
            arg_type = Int64
            default = 100000
        "--filename"
            arg_type = String
            default = ""
    end

    return parse_args(s)
end=#

function run!()
    # basic simulation parameters
    #   - radius: radius of particle
    #   - φ: packing fraction
    #   - f_propulsion: propulsion force (γ * f_propulsion is the propulsion speed)
    #   - D_rot: rotational diffusion constant (if using Active Brownian Particles)
    #       - α and τ are for Run-and-Tumble and Ornstein-Uhlenbeck
    #   - kT: temperature
    #   - num_steps: number of steps to run the simulation for
    #   - dt: timestep
    radius = 0.5

    φ = 0.6
    L = 120.0
    f_propulsion = 0.5
    D_rot = 0.001
    #α = 0.001
    #τ = 1000.0    
    kT = 0.0

    num_steps = 1000000
    dt = 0.01

    # information needed for saving simulation
    #   - filenames and directories for saving trajectories and movies (here we save the frames toward the end of and throughout simulation)
    #   - num_frames_to_save: number of frames to save
    #   - save_start: frame number to start recording simulation
    #   - save_interval: how often to record simulation
    traj_partial_filename = "TEST_OUTPUTS/trajectory_partial.out"
    animation_partial_filename = "TEST_OUTPUTS/animation_partial.mp4"
    traj_full_filename = "TEST_OUTPUTS/trajectory_full.out"
    animation_full_filename = "TEST_OUTPUTS/animation_full.mp4"
    num_frames_to_save = 300
    save_interval_partial = trunc(Int64, 0.5 / (f_propulsion * dt))
    save_start_partial = num_steps - save_interval_partial * num_frames_to_save
    save_interval_full = trunc(Int64, num_steps / num_frames_to_save)
    save_start_full = num_steps - save_interval_full * num_frames_to_save

    # initialize particles on a triangular lattice given a packing fraction and box size
    lattice_spacing = sqrt(2 * pi * radius^2 / (sqrt(3) * φ))
    Nx = trunc(Int64, L / lattice_spacing)
    initial_positions, box = initialize_triangular_lattice(lattice_spacing, (Nx, trunc(Int64, Nx / sqrt(3))))
    particles = Vector{Particle}()
    for position in initial_positions
        push!(particles, Particle(position, ActiveBrownian(f_propulsion, D_rot)))
        #push!(particles, Particle(position, RunAndTumble(f_propulsion, α)))
        #push!(particles, Particle(position, OrnsteinUhlenbeck(f_propulsion, τ)))
    end
    println("Number of particles: $(length(particles))")
    println("Packing fraction: $(φ)")
    println("Box: $(box)")

    # define cell lists for efficient neighbor finding
    cell_list, padding = CellList(particles, 2 * radius, box; approx_padding = 0.25)
    cell_list.update_interval = maximum([1, floor(Int64, 0.90 * padding / (f_propulsion * dt))])
    println("CellList update interval set to $(cell_list.update_interval)")
    
    # initialize repulsive interaction (for an elastic repulsion, though Lennard-Jones can also be used)
    #   - k: spring constant of the elastic repulsion
    k = 10.0
    hr = HarmonicRepulsion(k, 2 * radius, cell_list, box)

    # initialize integrator (Euler-Maruyama method)
    brownian = BrownianDynamics(particles, dt, kT, box)

    # collect all the pieces and initialize system for simulation
    system = System(particles, [hr], [cell_list], brownian)

    # for storing particle positions
    partial_trajectory = TrajectoryContainer(save_interval_partial, box; start = save_start_partial)
    full_trajectory = TrajectoryContainer(save_interval_full, box; start = save_start_full)

    # run and save simulation
    run_simulation!(system, [partial_trajectory, full_trajectory], num_steps)
    save!(partial_trajectory, traj_partial_filename)
    save!(full_trajectory, traj_full_filename)

    # functions for visualzation
    #   - visualize! pulls up a window with slider to view simulations
    #   - save_movie! saves simulation as a video
    #visualize!(partial_trajectory)
    #visualize!(full_trajectory)
    save_movie!(partial_trajectory, 30, animation_partial_filename)
    save_movie!(full_trajectory, 30, animation_full_filename)
end
run!()