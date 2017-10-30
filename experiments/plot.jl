using Gadfly
using Distributions
using Colors
using DataFrames
using Query

Gadfly.push_theme(Theme(major_label_font="Droid Sans",
                        minor_label_font="Droid Sans",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=14pt, key_position=:none,
                        default_color=colorant"#000000"))

colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a",
          colorant"#984ea3", colorant"#ff7f00", colorant"#ffff33"]
# colors = [colorant"#7fc97f", colorant"#beaed4", colorant"#fdc086",
#           colorant"#ffff99", colorant"#386cb0", colorant"#f0027f"]

function reducef(df, xmax)
    r = 1:length(df[:eval])
    if df[:eval][end] != xmax
        r = 1:(length(df[:eval])+1)
    end
    map(i->mapf(i, df, xmax), r)
end

function mapf(i::Int64, df, xmax::Int64)
    if i > length(df[:eval])
        return df[:fit][end] * ones(xmax - df[:eval][end])
    end
    lower = 0
    if i >= 2
        lower = df[:eval][i-1]
    end
    df[:fit][i] * ones(df[:eval][i] - lower)
end

function idf(ea::String, chromosome::String, mutation::String, crossover::String,
             distance::String)::Int64
    id = 1e4 * findfirst(["oneplus", "NEAT", "GA"] .== ea)
    id += 1e3 * findfirst(["CGP.CGPChromo", "CGP.EPCGPChromo", "CGP.RCGPChromo",
                           "CGP.PCGPChromo", "CGP.RPCGPChromo"] .== chromosome)
    id += 1e2 * findfirst(["CGP.mutate_genes", "CGP.mixed_node_mutate",
                           "CGP.mixed_subtree_mutate"] .== mutation)
    id += 1e1 * findfirst(["N/A", "CGP.single_point_crossover",
                           "CGP.proportional_crossover",
                           "CGP.random_node_crossover",
                           "CGP.aligned_node_crossover",
                           "CGP.output_graph_crossover",
                           "CGP.subgraph_crossover"] .== crossover)
    id += findfirst(["N/A", "CGP.functional_distance", "CGP.genetic_distance",
                     "CGP.positional_distance"] .== distance)
    id
end

function get_sweep_stats(log::String)
    res = readtable(log, header=false, separator=' ',
                    names=[:date, :time, :seed, :eval, :fit, :refit, :mean_fit,
                           :active, :nodes, :mean_nodes, :species, :ea, :chromosome,
                           :mutation, :crossover, :distance])

    xmax = maximum(res[:eval])

    filled = by(res, [:seed, :ea, :chromosome, :mutation, :crossover, :distance],
                df->reduce(vcat, reducef(df, xmax)))
    filled[:xs] = repeat(1:xmax, outer=Int64(size(filled,1)/xmax))
    stats = by(filled, [:xs, :ea, :chromosome, :mutation, :crossover, :distance],
               df->DataFrame(stds=std(df[:x1]), means=mean(df[:x1]), mins=minimum(df[:x1]),
                             maxs=maximum(df[:x1])))
    stats[:lower] = stats[:means]-0.5*stats[:stds]
    stats[:upper] = stats[:means]+0.5*stats[:stds]
    stats[:id] = idf.(stats[:ea], stats[:chromosome], stats[:mutation],
                      stats[:crossover], stats[:distance])
    stats
end

function get_smaller(stats::DataFrame)
    finals = by(stats, [:id], df->df[:means][end])
    sort!(finals, cols=[:x1], rev=true)
    ids = [finals[:id][1:3]; finals[:id][end]]
    if ~(11111 in ids)
        append!(ids, [11111])
    end

    x = @from i in stats begin
        @where i.id in ids
        @select {i.id, i.xs, i.means, i.stds, i.mins, i.maxs, i.lower, i.upper}
        @collect DataFrame
    end
    x
    i=isnan.(x[:stds]);
    x[:stds][i] = x[:means][i];
    x[:lower][i] = x[:means][i];
    x[:upper][i] = x[:means][i];
    x
end

function plot_sweep(stats::DataFrame, title::String, filename::String="training";
                    xmin=minimum(stats[:xs]),
                    xmax=maximum(stats[:xs]),
                    ymin=minimum(stats[:lower]),
                    ymax=maximum(stats[:upper]))
    plt = plot(stats, x="xs", y="means", ymin="lower", ymax="upper", color="id",
               Geom.line, Geom.ribbon,
               Scale.color_discrete_manual(colors...),
               Guide.title(title),
               Guide.xlabel("Evaluations"),
               Guide.ylabel("Fitness"),
               Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax))
    draw(PDF(string(filename, ".pdf"), 8inch, 6inch), plt)
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

function get_cmaes_results(log::String)
    res = readtable(log, header=false, separator=' ',
                    names=[:date, :time, :seed, :eval, :fit, :active, :nodes, :ea,
                           :chromo, :mutation, :crossover, :distance, :recurrency,
                           :input_mutation_rate, :output_mutation_rate, :node_mutation_rate,
                           :add_node_rate, :delete_node_rate, :add_mutation_rate,
                           :delete_mutation_rate, :speciation_thresh, :ga_elitism_rate,
                           :ga_crossover_rate, :ga_mutation_rate])
    res[:chromosome] = map(x->split(split(x, '.')[2], "Chromo")[1], res[:chromo])
    runs = by(res, [:ea, :chromosome, :mutation, :crossover, :distance, :recurrency,
                      :input_mutation_rate], df->maximum(df[:fit]))
    top_runs = by(runs, [:ea, :chromosome], df->sort!(df[:x1], rev=true)[1:20])

    res, top_runs
end

function plot_top_runs(top_runs::DataFrame, ti::String, filename::String;
                       ymin = floor(minimum(top_runs[:x1])*100)/100,
                       ymax = ceil(maximum(top_runs[:x1])*100)/100,
                       step = round((ymax-ymin)/5*10)/10)

    cgp = @from i in top_runs begin
        @where i.chromosome == "CGP"
        @select {i.x1, i.ea, i.chromosome}
        @collect DataFrame
    end
    pcgp = @from i in top_runs begin
        @where i.chromosome == "PCGP"
        @select {i.x1, i.ea, i.chromosome}
        @collect DataFrame
    end

    if step == 0
        step = (ymax-ymin)/5
    end

    yticks = [ymin:step:ymax;]

    plt1 = plot(layer(cgp, x=:ea, y=:x1, color=:ea, Geom.boxplot),
                Scale.color_discrete_manual(colors...),
                Coord.cartesian(ymin=ymin, ymax=ymax),
                Guide.yticks(ticks=yticks, label=true),
                Guide.ylabel("Fitness"), Guide.xlabel("CGP"));
    plt2 = plot(layer(pcgp, x=:ea, y=:x1, color=:ea, Geom.boxplot),
                Scale.color_discrete_manual(colors...),
                Coord.cartesian(ymin=ymin, ymax=ymax),
                Guide.yticks(ticks=yticks, label=false),
                Guide.ylabel(nothing), Guide.xlabel("PCGP"));

    plt = title(hstack(plt1, plt2), ti);
    draw(PDF(filename, 8inch, 6inch), plt)
    plt
end

