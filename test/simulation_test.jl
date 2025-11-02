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
    traj_filename = "TEST_OUTPUTS/trajectories.out"
    animation_filename = "TEST_OUTPUTS/animation.mp4"

    radius = 0.5

    phi = 0.5
    lattice_spacing = sqrt(2 * pi * radius^2 / (sqrt(3) * phi))
    Nx = 50
    initial_positions, box = initialize_triangular_lattice(lattice_spacing, (Nx, trunc(Int64, Nx / sqrt(3))))
    f_propulsion = 1.0
    D_rot = 0.001
    kT = 0.0

    num_steps = 100000
    save_interval = 25
    save_start = num_steps - save_interval * 600
    dt = 0.01

    k = 10.0

    particles = Vector{Particle}()
    for position in initial_positions
        push!(particles, Particle(position, ActiveBrownian(f_propulsion, D_rot)))
    end
    println("Number of particles: $(length(particles))")
    println("Packing fraction: $(phi)")

    cell_list, padding = CellList(particles, 2 * radius, box; approx_padding = 0.25)
    cell_list.update_interval = maximum([1, floor(Int64, 0.90 * padding / (f_propulsion * dt))])
    println("CellList update interval set to $(cell_list.update_interval)")
    lj = HarmonicRepulsion(k, 2 * radius, cell_list, box)

    brownian = BrownianDynamics(particles, dt, kT, box)

    system = System(particles, [lj], [cell_list], brownian)
    trajectories = Trajectories(save_interval, box; start = save_start)
    run_simulation!(system, trajectories, num_steps)
    save!(trajectories, traj_filename)

    #visualize!(trajectories)
    #save_movie!(trajectories, 30, animation_filename)
end
run!()