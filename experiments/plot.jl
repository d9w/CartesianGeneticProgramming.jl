using Gadfly
using Distributions
using Colors
using DataFrames

Gadfly.push_theme(Theme(major_label_font="Droid Sans",
                        minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=14pt, key_position=:right))

colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a",
          colorant"#984ea3", colorant"#ff7f00", colorant"#ffff33",
          colorant"#a65628", colorant"#f781bf"]

function get_raw_results(logdirs::Array{String}, xmax::Int64=25000, nruns::Int64=20)
    nlogs = length(logdirs)
    trains = DataFrame()
    xs = Array{Int64}(0)
    fits = Array{Float64}(0)
    logi = Array{Int64}(0)
    logs = Array{String}(0)

    for l in eachindex(logdirs)
        logdir = logdirs[l]
        lname = split(logdir, "/")[end]
        for i=0:(nruns-1)
            file = join([logdir, "/", string(i), ".log"])
            res = readdlm(file, skipstart=0, ' ')
            extra = 1
            if res[1, 3] != 0
                append!(xs, [0])
                append!(fits, res[1, 4])
                extra = 2
            end
            append!(xs, res[:, 3])
            append!(xs, [xmax])
            append!(fits, res[:, 4])
            append!(fits, [res[end, 4]])
            append!(logi, repeat([i], inner=size(res,1)+extra))
            append!(logs, repeat([lname], inner=size(res,1)+extra))
        end
    end

    trains[:xs] = xs
    trains[:fits] = fits
    trains[:logi] = logi
    trains[:logs] = logs
    trains[:mod_logi] = mod.(trains[:logi], 10)
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

function get_stats(trains::DataFrame)
    filled = by(trains, [:logi, :logs], df->reduce(
        vcat, map(i->df[:fits][i-1]*ones(df[:xs][i]-df[:xs][i-1]),
                  2:length(df[:xs]))));

    xmax = maximum(trains[:xs])
    filled[:xs] = repeat(1:xmax, outer=Int64(size(filled,1)/xmax))

    stats = by(filled, [:logs, :xs], df->DataFrame(
        stds=std(df[:x1]), means=mean(df[:x1]),
        mins=minimum(df[:x1]), maxs=maximum(df[:x1])))
    stats[:lower] = stats[:means]-stats[:stds]
    stats[:upper] = stats[:means]+stats[:stds]
    stats
end

function plot_points(trains::DataFrame, stats::DataFrame, title::String)
    plt = plot(layer(stats, x="xs", y="means", ymin="lower", ymax="upper",
                     color="logs", Geom.line, Geom.ribbon),
               layer(trains, x="xs", y="fits", color="logs", shape="mod_logi",
                     Geom.point),
               Scale.color_discrete_manual(colors...),
               Guide.title(title),
               Guide.xlabel("Evaluations"),
               Guide.ylabel("Fitness"))
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


