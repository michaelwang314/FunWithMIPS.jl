mutable struct CellList
    particles::Vector{Particle}
    start_index::Array{Int64, 2}
    next_index::Vector{Int64}

    cell_counts::SVector{2, Int64}
    cell_sizes::SVector{2, Float64}
    box::SVector{2, Float64}

    update_interval::Int64
    update_counter::Int64
end

function CellList(particles::Vector{Particle}, approx_cell_size::Float64, box::Vector{Float64}; 
                  update_interval::Int64 = 1, approx_padding::Float64 = 0.0)
    cell_counts = floor.(Int64, box ./ (approx_cell_size + approx_padding))
    cell_sizes = box ./ cell_counts

    start_index = -ones(Int64, cell_counts[1], cell_counts[2])
    next_index = -ones(Int64, length(particles))
    
    for (n, particle) in enumerate(particles)
        i, j = floor.(Int64, particle.position ./ cell_sizes) .+ 1

        if start_index[i, j] > 0
            next_index[n] = start_index[i, j]
        end
        start_index[i, j] = n
    end

    return CellList(particles, start_index, next_index, cell_counts, cell_sizes, box, update_interval, 0), minimum(cell_sizes .- approx_cell_size)
end

function update_neighbor_list!(cell_list::CellList)
    if (cell_list.update_counter += 1) == cell_list.update_interval
        fill!(cell_list.start_index, -1)

        @inbounds for n = 1 : length(cell_list.particles)
            x, y = cell_list.particles[n].position
            i = floor(Int64, x / cell_list.cell_sizes[1]) + 1
            j = floor(Int64, y / cell_list.cell_sizes[2]) + 1

            if cell_list.start_index[i, j] > 0
                cell_list.next_index[n] = cell_list.start_index[i, j]
            else
                cell_list.next_index[n] = -1
            end
            cell_list.start_index[i, j] = n
        end
        cell_list.update_counter = 0
    end
end