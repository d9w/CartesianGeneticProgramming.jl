using ArcadeLearningEnvironment
using Logging
using CGP

function play_atari(c::Chromosome, game::Game, id::String;
                    make_draw::Bool=false, folder::String="",
                    max_frames=18000)
    reset_game(game.ale)
    reward = 0.0
    frames = 0
    p_action = game.actions[1]
    while ~game_over(game.ale)
        output = process(c, get_rgb(game))
        action = game.actions[indmax(output)]
        reward += act(game.ale, action)
        if rand() < 0.25
            reward += act(game.ale, action) # repeat action here for seeding
        end
        if make_draw
            screen = draw(game)
            filename = string(folder, "/", @sprintf("frame_%06d.png", frames))
            save(filename, screen)
        end
        frames += 1
        if frames > max_frames
            Logging.debug(string("Termination due to frame count on ", id))
            break
        end
    end
    reward
end
