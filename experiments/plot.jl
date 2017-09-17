using Gadfly
using Distributions
using Colors
using DataFrames

Gadfly.push_theme(Theme(major_label_font="Droid Sans",
                        minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.3),
                        key_label_font_size=14pt, key_position=:right))

colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a",
          colorant"#984ea3", colorant"#ff7f00", colorant"#ffff33",
          colorant"#a65628", colorant"#f781bf"]

function get_julia_results(logdirs::Array{String}, xmax::Int64=50000, nruns::Int64=1)

    nlogs = length(logdirs)
    trains = zeros(nlogs, xmax, nruns)

    for l in eachindex(logdirs)
        logdir = logdirs[l]
        for i=0:(nruns-1)
            file = join([logdir, "/", string(i), ".log"])
            res = readdlm(file, skipstart=1, ' ')[:, 3:end]
            cmin = 1
            for c=1:size(res)[1]
                cur = Int(res[c,1])
                if cur > xmax
                    break
                end
                trains[l, cmin:cur, i+1] = res[c,2]
                cmin = cur + 1
            end
            trains[l, cmin:xmax, i+1] = res[size(res)[1],2]
        end
    end

    trains
end

function get_cpp_results(logdirs::Array{String}, xmax::Int64=50000, nruns::Int64=20)

    nlogs = length(logdirs)
    trains = zeros(nlogs, xmax, nruns)
    tests = zeros(nlogs, xmax, nruns)
    sizes = zeros(nlogs, xmax, nruns)

    for l in eachindex(logdirs)
        logdir = logdirs[l]

        for i=0:(nruns-1)
            file = join([logdir, "/", string(i), ".log"])
            res = readcsv(file)
            if res[1,1] != 1.0
                println(file)
                error(file)
            end
            cmin = 1
            for c=1:size(res)[1]
                cur = Int(res[c,1])
                if cur > xmax
                    break
                end
                trains[l, cmin:cur, i+1] = res[c,2]
                tests[l, cmin:cur, i+1] = res[c,3]
                sizes[l, cmin:cur, i+1] = res[c,4]
                cmin = cur + 1
            end
            trains[l, cmin:xmax, i+1] = res[size(res)[1],2]
            tests[l, cmin:xmax, i+1] = res[size(res)[1],3]
            sizes[l, cmin:xmax, i+1] = res[size(res)[1],4]
        end
    end

    trains, tests, sizes
end

function plot_training(trains::Array{Float64}, valids::Array{Float64},
                      labels::Array{String}, title::String)

    nlogs, xmax, nruns = size(trains)
    layers = Array{Array{Gadfly.Layer,1}}(nlogs*2)

    for l in 1:nlogs
        mtrain = mean(trains[l,:,:], 2)
        strain = std(trains[l,:,:], 2)
        mvalid = mean(valids[l,:,:], 2)
        svalid = std(valids[l,:,:], 2)

        layers[l*2-1] = layer(x=1:xmax, y=mtrain,
                              ymin=mtrain-0.5*strain, ymax=mtrain+0.5*strain,
                              Geom.line, Geom.ribbon,
                              style(default_color=colors[l], line_style=:solid))
        layers[l*2] = layer(x=1:xmax, y=mvalid,
                            ymin=mvalid-0.5*svalid, ymax=mvalid+0.5*svalid,
                            Geom.line, Geom.ribbon,
                            style(default_color=colors[l], line_style=:dash,
                                  lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.1)))
    end

    plt = plot(layers..., Guide.title(title),
               Guide.xlabel("Generation"), Guide.ylabel("Fitness"),
               # Guide.yticks(ticks=[0:0.2:1.0;]),
               # Scale.x_continuous(labels=x -> @sprintf("%de3", x/1e3)),
               Guide.manual_color_key(
                   "Legend", labels, colors[1:nlogs]));
    draw(PDF("training.pdf", 8inch, 6inch), plt)
    plt
end

function plot_tests(tests::Array{Float64}, labels::Array{String}, title::String)

    d = DataFrame()
    d[:Accuracy] = tests'[:]
    d[:Dataset] = repeat(labels, inner=20)
    d[:Mean] = repeat(mean(tests,2)[:], inner=20)

    plt = plot(
        # layer(d, x="Dataset", y="Mean", Geom.point,
                     # style(default_color=colorant"#000000")),
               layer(d, x="Dataset", y="Accuracy", Geom.violin,
                     style(default_color=colorant"#000000")),
               Guide.yticks(ticks=[0:0.2:1.0;]),
               Guide.ylabel("Accuracy"), Guide.xlabel(nothing),
               Guide.title(title));
    draw(PDF("tests.pdf", 8inch, 6inch), plt)
    plt
end


