using ArcadeLearningEnvironment
using CGP
using Logging
# using Images

function play_atari(c::Chromosome, game::Game)
    reset_game(game.ale)
    reward = 0
    # frames = 0
    while ~game_over(game.ale)
        act(game.ale, game.actions[2])
        output = process(c, get_inputs(game))
        reward += act(game.ale, game.actions[indmax(output)])
        # frames += 1
        # save(@sprintf("breakout/frame_%05d.png", frames), draw(game))
    end
    reward
end

seed = 0
log = "log"
ea = oneplus
ctype = CGPChromo
if length(ARGS) > 0; seed = parse(Int64, ARGS[1]); end
if length(ARGS) > 1; log = ARGS[2]; end
if length(ARGS) > 2; ea = eval(parse(ARGS[3])); end
if length(ARGS) > 3; ctype = eval(parse(ARGS[4])); end

CGP.Config.init("cfg/base.yaml")
if ctype == MTPCGPChromo
    CGP.Config.init("cfg/mtcgp.yaml")
else
    CGP.Config.init("cfg/classic.yaml")
end

Logging.configure(filename=log, level=INFO)
Logging.info("I: $seed $ea $ctype")
srand(seed)

game = Game("breakout")
nin = length(get_inputs(game))
nout = length(game.actions)
fit = x->play_atari(x, game)
maxfit, best = ea(ctype, nin, nout, fit)
close!(game)
