function visualize!(trajectories::Trajectories)
    max_time = length(trajectories)
    num_particles = length(trajectories.history[1])

    positions = [[Point2f(particle.position)] for particle in trajectories.history[1]]
    for t = 2 : max_time
        for (p, particle) in enumerate(trajectories.history[t])
            push!(positions[p], Point2f(particle.position))
        end
    end

    window = Figure(size = (800, 600))
    scene_axis = Axis(window[2, 1])
    hidedecorations!(scene_axis)
    limits!(scene_axis, 0, trajectories.box[1], 0, trajectories.box[2])

    controls_grid = GridLayout(window[1, 1], tell_width = false)
    time_slider = Slider(controls_grid[1, 1], range = 1 : 1 : max_time, startvalue = 1, update_while_dragging = true)
    Label(controls_grid[1, 2], "Play")
    speed_menu = Menu(controls_grid[1, 3], options = ["x1", "x2", "x4", "x8", "x16"], default = "x1", width = 50)
    play_button = Toggle(controls_grid[1, 4], active = false)

    speed = Observable{Int64}(1)
    on(speed_menu.selection) do s
        speed[] = parse(Int64, chop(s, head = 1, tail = 0))
    end

    on(fig.scene.events.tick) do tick
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
        object = lift(time_slider.value) do t
            Circle(positions[p][t], 0.5f0)
        end
        poly!(scene_axis, object, color = :blue)
    end

    wait(display(window))
end