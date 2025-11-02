module FunWithMIPS
    using StaticArrays
    using DataStructures
    using Serialization
    using GLMakie

    export Propulsion, ActiveBrownian, RunAndTumble, OrnsteinUhlenbeck, Particle
    export CellList, update_neighbor_list!
    export Interaction, LennardJones, HarmonicRepulsion, compute_forces!, correct_for_periodicity
    export Integrator, BrownianDynamics, update_particles!, get_and_update_active_force
    export System, TrajectoryContainer, initialize_triangular_lattice, run_simulation!, save!, load
    export visualize!, save_movie!

    """
    Flush output so that jobs can be monitored on a cluster.
    """
    @inline println(args...) = println(stdout, args...)
    @inline function println(io::IO, args...)
        Base.println(io, args...)
        flush(io)
    end

    include("particles.jl")
    include("neighbor_lists.jl")
    include("interactions.jl")
    include("integrators.jl")
    include("system.jl")

    include("visualizer/render.jl")
end
