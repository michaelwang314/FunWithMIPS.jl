module FunWithMIPS
    using StaticArrays
    using DataStructures
    using GLMakie

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
end
