struct System
    particles::Vector{Particle}
    interactions::Vector{<:Interaction}
    neighbor_lists::Vector{CellList}
    integrator::Integrator
end

mutable struct Trajectories
    history::Vector{Vector{Particle}}

    start::Int64
    period::Int64

    box::SVector{2, Float64}
end
Trajectories(period::Int64, box::Vector{Float64}; start::Int64 = 1) = Trajectories(Vector{Vector{Particle}}(), start, period, box)

function run_simulation!(system::System, trajectories::Union{Trajectories, Nothing}, num_steps::Int64; message_interval::Union{Float64, Nothing} = 10.0)
    prev_step = 0
    time_elapsed = 0.0
    interval_start = time()
    for step = 1 : num_steps
        for interaction in system.interactions
            compute_forces!(interaction)
        end

        update_particles!(system.integrator)
        
        for neighbor_list in system.neighbor_lists
            update_neighbor_list!(neighbor_list)
        end

        if !isnothing(trajectories) && (step - trajectories.start) % trajectories.period == 0
            push!(trajectories.history, deepcopy(system.particles))
        end

        interval_time = time() - interval_start
        if !isnothing(message_interval) && (interval_time > message_interval || step == num_steps)
            time_elapsed += interval_time
            rate = (step - prev_step) / interval_time
            println(hr_min_sec(time_elapsed), " | ",
                    step, "/", num_steps, " (", round(step / num_steps * 100, digits = 1), "%) | ",
                    round(rate, digits = 1), " steps/s | ",
                    hr_min_sec((num_steps - step) / rate))
            prev_step = step
            interval_start = time()
        end
    end
    println("Average steps/s: ", round(num_steps / time_elapsed, digits = 1))
end
run_simulation!(system::System, num_steps::Int64; message_interval::Union{Float64, Nothing} = 10.0) = run_simulation!(system, save_to, num_steps; message_interval = message_interval)

function initialize_triangular_lattice(lattice_spacing::Float64, dims::Tuple{Int64, Int64})
    Nx, Ny = dims
    
    positions = Vector{Vector{Float64}}()
    unit_cell = [[0.0, 0.0], [lattice_spacing / 2, sqrt(3) * lattice_spacing / 2]]
    for nx = 0 : Nx - 1, ny = 0 : Ny - 1
        push!(positions, unit_cell[1] + [lattice_spacing * nx, sqrt(3) * lattice_spacing * ny])
        push!(positions, unit_cell[2] + [lattice_spacing * nx, sqrt(3) * lattice_spacing * ny])
    end

    return positions, [lattice_spacing * Nx, sqrt(3) * lattice_spacing * Ny]
end

@inline function hr_min_sec(time::Float64)
    hours = trunc(Int64, time / 3600.0)
    minutes = trunc(Int64, mod(time, 3600.0) / 60.0)
    seconds = trunc(Int64, mod(time, 60.0))

    return string(hours < 10 ? "0" : "", hours, 
                  minutes < 10 ? ":0" : ":", minutes, 
                  seconds < 10 ? ":0" : ":", seconds)
end

function save!(trajectories::Trajectories, filename::String)
    if !isdir(dirname(filename))
        mkpath(dirname(filename))
    end

    open(filename, "w") do io
        serialize(io, trajectories)
    end
    println("$(filename) saved!")
end

function load(filename::String)
    object = begin
        open(filename, "r") do io
            deserialize(io)
        end
    end
    println("$(filename) loaded!")

    return object
end