using FunWithMIPS
#=using ArgParse

function parse_cmd()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--"
            arg_type = Float64
            default = 1.0
        "--numlinkers"
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

    initial_positions, box = initialize_triangular_lattice(1.5, (20, 12))
    f_propulsion = 1.0
    D_rot = 0.01
    kT = 0.0

    num_steps = 100000
    save_every = 100
    dt = 0.01

    k = 10.0
    σ = 0.5

    particles = Vector{Particle}()
    for position in initial_positions
        push!(particles, Particle(position, ActiveBrownian(f_propulsion, D_rot)))
    end

    cell_list, padding = CellList(particles, 2 * σ, box)
    cell_list.update_interval = maximum(1, floor(Int64, 0.90 * padding / (f_propulsion * dt)))
    println("CellList update interaval set to $(cell_list.update_interval)")
    lj = HarmonicRepulsion(k, 2 * σ, cell_list, box)

    brownian = BrownianDynamics(particles, dt, kT, box)

    system = System(particles, [lj], [cell_list], brownian)
    trajectories = Trajectories(save_every, box)
    run_simulation!(system, trajectories, num_steps)
    save!(trajectories, traj_filename)

    visualize!(trajectories)
end
run!()