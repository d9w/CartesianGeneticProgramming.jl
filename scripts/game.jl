# The following was inspired by AtariAlgos.jl
# Modifications by Dennis Wilson @d9w

using Colors
using ImageCore
using ImageTransformations
export
    Game,
    close!,
    draw,
    get_inputs,
    get_rgb

struct Game
    ale::ALEPtr
    width::Int
    height::Int
    actions::Array{Int32}
end

function Game(romfile::String, seed::Int)
    ale = ALE_new()
    setInt(ale, "random_seed", Cint(seed))
    loadROM(ale, romfile)
    w = getScreenWidth(ale)
    h = getScreenHeight(ale)
    actions = getMinimalActionSet(ale)
    Game(ale, w, h, actions)
end

function close!(game::Game)
    ALE_del(game.ale)
end

function draw(game::Game)
    rawscreen = getScreenRGB(game.ale)
    colorview(RGB, Float64.(reshape(rawscreen/256.,
                                    (3, game.width, game.height))))';
end

function get_inputs(game::Game)
    screen = getScreen(game.ale)/(0xff*1.0)
    screen = reshape(screen, (game.width, game.height))'
    # imresize(screen, (42, 32))/256.
    screen
end

function get_rgb(game::Game)
    # TODO: speed this up
    # rawscreen = Array{Cuchar}(undef, game.width * game.height * 3)
    rawscreen = getScreenRGB(game.ale)
    rgb = Float64.(reshape(rawscreen/256., (3, game.width, game.height)));
    [Array{Float64}(rgb[i,:,:]) for i in 1:3]
end

function get_ram(game::Game)
    getRAM(game.ale) ./ typemax(UInt8)
end
