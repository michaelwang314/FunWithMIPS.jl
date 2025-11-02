abstract type Integrator end

struct BrownianDynamics <: Integrator
    particles::Vector{Particle}

    dt::Float64
    kT::Float64

    box::SVector{2, Float64}
end

function update_particles!(integrator::BrownianDynamics)
    Threads.@threads for particle in integrator.particles
        fa_x, fa_y = get_and_update_active_force(particle.propulsion, integrator.dt)

        dt_scaled = integrator.dt / particle.γ_trans
        if integrator.kT > 0.0
            thermal_amp = sqrt(2 * integrator.kT * dt_scaled)
            particle.position[1] += dt_scaled * (particle.force[1] + fa_x) + thermal_amp * randn()
            particle.position[2] += dt_scaled * (particle.force[2] + fa_y) + thermal_amp * randn()
        else
            particle.position[1] += dt_scaled * (particle.force[1] + fa_x)
            particle.position[2] += dt_scaled * (particle.force[2] + fa_y)
        end

        fill!(particle.force, 0.0)

        if !(0.0 <= particle.position[1] < integrator.box[1])
            particle.position[1] = mod(particle.position[1], integrator.box[1])
        end
        if !(0.0 <= particle.position[2] < integrator.box[2])
            particle.position[2] = mod(particle.position[2], integrator.box[2])
        end
    end
end

@inline function get_and_update_active_force(propulsion::ActiveBrownian, dt::Float64)
    fx, fy = propulsion.force
    
    if propulsion.D_rot > 0.0
        δθ = sqrt(2 * propulsion.D_rot * dt) * randn()
        sin, cos = sincos(δθ)

        propulsion.force[1] = fx * cos - fy * sin
        propulsion.force[2] = fx * sin + fy * cos
    end

    return fx, fy
end

@inline function get_and_update_active_force(propulsion::RunAndTumble, dt::Float64)
    fx, fy = propulsion.force

    propulsion.τ_run -= dt
    if propulsion.τ_run < 0.0
        propulsion.τ_run = -log(rand()) / propulsion.α
        δθ = 2 * pi * rand()
        sin, cos = sincos(δθ)

        propulsion.force[1] = fx * cos - fy * sin
        propulsion.force[2] = fx * sin + fy * cos
    end

    return fx, fy
end

@inline function get_and_update_active_force(propulsion::OrnsteinUhlenbeck, dt::Float64)
    fx, fy = propulsion.force

    dt_scaled = dt / propulsion.τ
    noise_amp = sqrt(propulsion.strength * dt_scaled)
    propulsion.force[1] = fx - dt_scaled * fx + noise_amp * randn()
    propulsion.force[2] = fy - dt_scaled * fy + noise_amp * randn()

    return fx, fy
end

