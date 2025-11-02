abstract type Propulsion end

mutable struct ActiveBrownian <: Propulsion
    force::MVector{2, Float64}

    D_rot::Float64
end
function ActiveBrownian(strength::Float64, D_rot::Float64)
    θ = 2 * pi * rand()
    force = strength * [cos(θ), sin(θ)]

    return ActiveBrownian(force, D_rot)
end

mutable struct RunAndTumble <: Propulsion
    force::MVector{2, Float64}
    τ_run::Float64

    α::Float64
end
function RunAndTumble(strength::Float64, α::Float64)
    θ = 2 * pi * rand()
    τ_run = -log(rand()) / α
    force = strength * [cos(θ), sin(θ)]

    return RunAndTumble(force, τ_run, α)
end

mutable struct OrnsteinUhlenbeck <: Propulsion
    force::MVector{2, Float64}
    
    strength::Float64
    τ::Float64
end
function OrnsteinUhlenbeck(strength::Float64, τ::Float64)
    force = strength / sqrt(2) * randn(2)

    return OrnsteinUhlenbeck(force, strength, τ)
end

struct Particle
    position::MVector{2, Float64}
    force::MVector{2, Float64}

    propulsion::Union{<:Propulsion, Nothing}

    γ_trans::Float64
end
Particle(position::Vector{Float64}, propulsion::Union{<:Propulsion, Nothing}) = Particle(position, [0.0, 0.0], propulsion, 1.0)