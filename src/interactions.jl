abstract type Interaction end

struct LennardJones <: Interaction
    ϵ::Float64
    σ::Float64
    r_cut::Float64

    neighbor_list::CellList
    box::SVector{2, Float64}
end

struct HarmonicRepulsion <: Interaction
    k::Float64
    r_cut::Float64

    neighbor_list::CellList
    box::SVector{2, Float64}
end

function compute_forces!(lj::LennardJones)
    Threads.@threads for particle in lj.neighbor_list.particles
        x, y = particle.position
        i = floor(Int64, x / lj.neighbor_list.cell_sizes[1])
        j = floor(Int64, y / lj.neighbor_list.cell_sizes[2])

        @inbounds for Δi = -1 : 1, Δj = -1 : 1
            iΔi = mod(i + Δi, lj.neighbor_list.cell_counts[1]) + 1
            jΔj = mod(j + Δj, lj.neighbor_list.cell_counts[2]) + 1

            index = lj.neighbor_list.start_index[iΔi, jΔj]
            while index > 0
                neighbor = lj.neighbor_list.particles[index]
                Δx, Δy = correct_for_periodicity(x - neighbor.position[1], y - neighbor.position[2], lj.box)
                
                Δr² = Δx^2 + Δy^2
                if 0.0 < Δr² < lj.r_cut^2
                    pow = (lj.σ / Δr²)^3
                    coef = lj.ϵ * (48.0 * pow - 24.0) * pow / Δr²

                    particle.force[1] += coef * Δx
                    particle.force[2] += coef * Δy
                end
                index = lj.neighbor_list.next_index[index]
            end
        end
    end
end

function compute_forces!(hr::HarmonicRepulsion)
    Threads.@threads for particle in hr.neighbor_list.particles
        x, y = particle.position
        i = floor(Int64, x / hr.neighbor_list.cell_sizes[1])
        j = floor(Int64, y / hr.neighbor_list.cell_sizes[2])

        @inbounds for Δi = -1 : 1, Δj = -1 : 1
            iΔi = mod(i + Δi, hr.neighbor_list.cell_counts[1]) + 1
            jΔj = mod(j + Δj, hr.neighbor_list.cell_counts[2]) + 1

            index = hr.neighbor_list.start_index[iΔi, jΔj]
            while index > 0
                neighbor = hr.neighbor_list.particles[index]
                Δx, Δy = correct_for_periodicity(x - neighbor.position[1], y - neighbor.position[2], hr.box)
                
                Δr² = Δx^2 + Δy^2
                if 0.0 < Δr² < hr.r_cut^2
                    Δr = sqrt(Δr²)
                    coef = hr.k * (hr.r_cut / Δr - 1)

                    particle.force[1] += coef * Δx
                    particle.force[2] += coef * Δy
                end
                index = hr.neighbor_list.next_index[index]
            end
        end
    end
end

@inline function correct_for_periodicity(Δx::Float64, Δy::Float64, box::SVector{2, Float64})
    if abs(Δx) > 0.5 * box[1]
        Δx -= sign(Δx) * box[1]
    end
    if abs(Δy) > 0.5 * box[2]
        Δy -= sign(Δy) * box[2]
    end

    return Δx, Δy
end