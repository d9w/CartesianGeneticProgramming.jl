using ArcadeLearningEnvironment
using CGP
using Logging
using Images


function play_qbert()
    game = Game("qbert", 0)
    reward = 0
    frames = 0
    while ~game_over(game.ale)
        inputs = get_inputs(game)
        action = 0
        if inputs[13076] > 0
            action = 2
        end
        reward += act(game.ale, Cint(action))
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

function random_breakout_expert(nin::Int64, nout::Int64)
    # TODO: debug / rewrite without segmentation
    c = HPCGPChromo(nin, nout)
    pos = sort!(rand(15))
    c.genes[1] = pos[1]
    for out in nin+(1:nout)
        c.genes[out] = pos[15]
    end
    c.genes[nin+4] = pos[13]
    c.genes[nin+5] = pos[14]
    CGP.set_genes!(c, nin+1, [
        pos[2], rand(), rand(), CGP.Config.func2f(CGP.Config.f_const), 0.6])
    CGP.set_genes!(c, nin+2, [
        pos[3], rand(), rand(), CGP.Config.func2f(CGP.Config.f_const), 0.8])
    CGP.set_genes!(c, nin+3, [
        pos[4], pos[1], rand(), CGP.Config.func2f(CGP.Config.f_felzenszwalb), 0.8])
    CGP.set_genes!(c, nin+4, [
        pos[5], pos[4], pos[2], CGP.Config.func2f(CGP.Config.f_gt), rand()])
    CGP.set_genes!(c, nin+5, [
        pos[6], pos[4], pos[3], CGP.Config.func2f(CGP.Config.f_lt), rand()])
    CGP.set_genes!(c, nin+6, [
        pos[7], pos[4], pos[3], CGP.Config.func2f(CGP.Config.f_gt), rand()])
    CGP.set_genes!(c, nin+7, [
        pos[8], pos[5], pos[6], CGP.Config.func2f(CGP.Config.f_and), rand()])
    CGP.set_genes!(c, nin+8, [
        pos[9], pos[7], rand(), CGP.Config.func2f(CGP.Config.f_com), rand()])
    CGP.set_genes!(c, nin+9, [
        pos[10], pos[8], rand(), CGP.Config.func2f(CGP.Config.f_com), rand()])
    CGP.set_genes!(c, nin+10, [
        pos[11], pos[9], rand(), CGP.Config.func2f(CGP.Config.f_last), rand()])
    CGP.set_genes!(c, nin+11, [
        pos[12], pos[10], rand(), CGP.Config.func2f(CGP.Config.f_last), rand()])
    CGP.set_genes!(c, nin+12, [
        pos[13], pos[11], pos[12], CGP.Config.func2f(CGP.Config.f_gt), rand()])
    CGP.set_genes!(c, nin+13, [
        pos[14], pos[11], pos[12], CGP.Config.func2f(CGP.Config.f_lt), rand()])
    CGP.set_genes!(c, nin+14, [
        pos[15], rand(), rand(), CGP.Config.func2f(CGP.Config.f_zeros), rand()])
    HPCGPChromo(c.genes, nin, nout)
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
            action = 5
        else
            action = 4
        end
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
