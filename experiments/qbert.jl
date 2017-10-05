using ArcadeLearningEnvironment
using CGP
using Logging
using Images

function play_qbert()
    game = Game("qbert")
    reward = 0
    frames = 0
    while ~game_over(game.ale)
        inputs = get_inputs(game)
        action = 0
        if inputs[13076] > 0
            action = 2
        end
        reward += act(game.ale, action)
        frames += 1
        screen = draw(game)
        for pixel_i in -1:1:1
            for pixel_j in -1:1:1
                if ~(pixel_i == 0 && pixel_j == 0)
                  screen[56+pixel_i, 63+pixel_j] = RGB{Float64}(1.0, 1.0, 1.0)
                end
            end
        end
        save(@sprintf("qbert/frame_%06d.png", frames), screen)
    end
    close!(game)
    reward
end
