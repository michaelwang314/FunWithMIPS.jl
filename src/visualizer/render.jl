function visualize!(trajectories::Trajectories)
    max_time = length(trajectories.history)
    num_particles = length(trajectories.history[1])

    positions = [[Point2f(particle.position)] for particle in trajectories.history[1]]
    for t = 2 : max_time
        for (p, particle) in enumerate(trajectories.history[t])
            push!(positions[p], Point2f(particle.position))
        end
    end

    time = Observable{Int64}(1)

    window = Figure(size = (800, 600))
    scene_axis = Axis(window[2, 1], aspect = trajectories.box[1] / trajectories.box[2], title = @lift("t = $($time)"))
    hidedecorations!(scene_axis)
    limits!(scene_axis, 0, trajectories.box[1], 0, trajectories.box[2])

    controls_grid = GridLayout(window[1, 1], tell_width = false)
    time_slider = Slider(controls_grid[1, 1], range = 1 : 1 : max_time, startvalue = 1, update_while_dragging = true)
    Label(controls_grid[1, 2], "Play")
    speed_menu = Menu(controls_grid[1, 3], options = ["x1", "x2", "x4", "x8", "x16"], default = "x1", width = 50)
    play_button = Toggle(controls_grid[1, 4], active = false)

    on(time_slider.value) do t
        time[] = t
    end

    speed = Observable{Int64}(1)
    on(speed_menu.selection) do s
        speed[] = parse(Int64, chop(s, head = 1, tail = 0))
    end
    on(window.scene.events.tick) do tick
        play_button.active[] || return
        if time_slider.value[] + speed[] < max_time
            time_slider.value[] += speed[]
            set_close_to!(time_slider, time_slider.value[])
        else
            time_slider.value[] = 1
            set_close_to!(time_slider, 1)
        end
    end

    for p = 1 : num_particles
        poly!(scene_axis, @lift(Circle(positions[p][$time], 0.5f0)), color = (:blue, 0.5), strokecolor = :black, strokewidth = 1.5)
    end

    wait(display(window))
end

function save_movie!(trajectories::Trajectories, framerate::Int64, filename::String)
    max_time = length(trajectories.history)
    num_particles = length(trajectories.history[1])

    positions = [[Point2f(particle.position)] for particle in trajectories.history[1]]
    for t = 2 : max_time
        for (p, particle) in enumerate(trajectories.history[t])
            push!(positions[p], Point2f(particle.position))
        end
    end

    time = Observable{Int64}(1)

    window = Figure(size = (800, 600))
    scene_axis = Axis(window[1, 1], aspect = trajectories.box[1] / trajectories.box[2], title = @lift("t = $($time)"))
    hidedecorations!(scene_axis)
    limits!(scene_axis, 0, trajectories.box[1], 0, trajectories.box[2])

    for p = 1 : num_particles
        poly!(scene_axis, @lift(Circle(positions[p][$time], 0.5f0)), color = (:blue, 0.5), strokecolor = :black, strokewidth = 1.5)
    end

    println("Saving movie to $(filename)")
    record(window, filename, 1 : 1 : max_time; framerate = framerate) do t
        time[] = t
        println("$(t)/$(max_time) done")
    end
end