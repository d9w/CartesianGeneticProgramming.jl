using ArcadeLearningEnvironment
# using CGP
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

function f_indmax(x::Array)
    i = indmax(x)
    if ndims(i) > 1
        i = ind2sub(x, i)
    end
    i ./ size(x)
end

f_indmax(x::Array, y::Array) = f_indmax(x)

function f_com(x::Array)
    # TODO: get coordinates of center of mass of array
end

function get_breakout_action(inputs::Array{Float64})
    action = 2
    bullet_box = inputs[94:189, 9:152];
    slider_box = inputs[190:192, 9:152];
    bullet_poses = find(bullet_box)
    slider_poses = find(slider_box)
    if length(bullet_poses) > 0 && length(slider_poses) > 0
        bullet_pos = mean(map(x->ind2sub(bullet_box, x)[2], bullet_poses))
        slider_pos = mean(map(x->ind2sub(slider_box, x)[2], slider_poses))
        if bullet_pos < slider_pos
            action = game.actions[5]
        else
            action = game.actions[4]
        end
    end
    if f_indmax(bullet_box) < f_indmax(slider_box)
        action = 5
    else
        action = 4
    end
    action
end

function play_breakout()
    game = Game("breakout")
    reward = 0
    frames = 0
    act(game.ale, game.actions[2])
    life = lives(game.ale)
    while ~game_over(game.ale)
        action = game.actions[2]
        if lives(game.ale) < life
            life = lives(game.ale)
        else
            inputs = get_inputs(game)
            action = game.actions[get_breakout_action(inputs)]
        end
        reward += act(game.ale, action)
        frames += 1
        screen = draw(game)
        save(@sprintf("breakout/frame_%06d.png", frames), screen)
    end
    close!(game)
    reward
end


# ACTION_MEANING = {
#     0 : "NOOP",
#     1 : "FIRE",
#     2 : "UP",
#     3 : "RIGHT",
#     4 : "LEFT",
#     5 : "DOWN",
#     6 : "UPRIGHT",
#     7 : "UPLEFT",
#     8 : "DOWNRIGHT",
#     9 : "DOWNLEFT",
#     10 : "UPFIRE",
#     11 : "RIGHTFIRE",
#     12 : "LEFTFIRE",
#     13 : "DOWNFIRE",
#     14 : "UPRIGHTFIRE",
#     15 : "UPLEFTFIRE",
#     16 : "DOWNRIGHTFIRE",
#     17 : "DOWNLEFTFIRE",
# }
