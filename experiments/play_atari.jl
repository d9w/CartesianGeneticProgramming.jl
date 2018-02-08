using ArcadeLearningEnvironment
using CGP

function play_atari(c::Chromosome, game::Game, id::String;
                    make_draw::Bool=false, folder::String="",
                    max_frames=18000, max_act_count=1000)
    reset_game(game.ale)
    reward = 0.0
    frames = 0
    p_action = game.actions[1]
    act_count = 0
    while ~game_over(game.ale)
        output = process(c, get_rgb(game))
        action = game.actions[indmax(output)]
        if action == p_action
            act_count += 1
        else
            p_action = action
            act_count = 0
        end
        if act_count > max_act_count
            println("Termination due to repetitive action ", id)
            return -Inf
        end
        reward += act(game.ale, action)
        frames += 1
        if rand() < 0.25
            reward += act(game.ale, action)
        end
        if make_draw
            screen = draw(game)
            filename = string(folder, "/", @sprintf("frame_%06d.png", frames))
            save(filename, screen)
        end
        if frames > max_frames
            println("Termination due to frame count on ", id)
            break
        end
    end
    if make_draw
        draw_graph(c, string(folder, "/graph.pdf"))
    end
    reward
end
