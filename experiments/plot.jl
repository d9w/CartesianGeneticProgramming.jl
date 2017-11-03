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

all_params = [:recurrency, :input_mutation_rate, :output_mutation_rate,
              :node_mutation_rate, :add_node_rate, :delete_node_rate,
              :add_mutation_rate, :delete_mutation_rate,
              :speciation_thresh, :ga_elitism_rate, :ga_crossover_rate,
              :ga_mutation_rate, :f_mutation, :f_crossover, :f_distance]

ea_params = [:recurrency, :input_mutation_rate, :output_mutation_rate,
             :node_mutation_rate, :add_node_rate, :delete_node_rate,
             :add_mutation_rate, :delete_mutation_rate,
             :f_mutation]

ga_params = [:recurrency, :input_mutation_rate, :output_mutation_rate,
             :node_mutation_rate, :add_node_rate, :delete_node_rate,
             :add_mutation_rate, :delete_mutation_rate, :ga_elitism_rate,
             :ga_crossover_rate, :ga_mutation_rate, :f_mutation, :f_crossover]

neat_params = [:recurrency, :input_mutation_rate, :output_mutation_rate,
               :node_mutation_rate, :add_node_rate, :delete_node_rate,
               :add_mutation_rate, :delete_mutation_rate,
               :speciation_thresh, :ga_crossover_rate,
               :ga_mutation_rate, :f_mutation, :f_crossover, :f_distance]

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
    maxeval = maximum(res[:eval])
    runs = @from i in res begin
        @where i.eval == maxeval
        @select i
        @collect DataFrame
    end
    mutations = sort!(unique(runs[:mutation]))
    println(mutations)
    runs[:f_mutation] = indexin(runs[:mutation], mutations)./(1.0*length(mutations))
    crossovers = sort!(unique(runs[:crossover]))
    println(crossovers)
    runs[:f_crossover] = indexin(runs[:crossover], crossovers)./(1.0*length(crossovers))
    distances = sort!(unique(runs[:distance]))
    println(distances)
    runs[:f_distance] = indexin(runs[:distance], distances)./(1.0*length(distances))
    eas = unique(runs[:ea])
    chromosomes = unique(runs[:chromosome])

    top_runs = DataFrame()
    for ea in eas
        for chromo in chromosomes
            group = @from i in runs begin
                @where i.ea == ea && i.chromosome == chromo
                @select i
                @collect DataFrame
            end
            sort!(group, cols=[:fit], rev=true)
            group[:id] = 1:size(group, 1)
            top_group = @from i in group begin
                @where i.id <= 20
                @select i
                @collect DataFrame

            end
            top_runs = vcat(top_runs, top_group)
        end
    end

    runs, top_runs
end

function plot_top_runs(top_runs::DataFrame, ti::String, filename::String;
                       ymin = floor(minimum(top_runs[:fit])*100)/100,
                       ymax = ceil(maximum(top_runs[:fit])*100)/100,
                       step = round((ymax-ymin)/5*10)/10)

    yticks = [ymin:step:ymax;]

    plt = plot(top_runs, xgroup=:chromosome,
               x=:ea, y=:fit, color=:ea,
               Geom.subplot_grid(Coord.cartesian(ymin=ymin, ymax=ymax),
                                 Geom.boxplot),
               Scale.color_discrete_manual(colors...),
               Guide.ylabel(nothing), Guide.xlabel(nothing),
               Guide.title(ti))
    draw(PDF(string(filename, "_fitness.pdf"), 8inch, 6inch), plt)
    plt
end

function plot_top_params(top_runs::DataFrame, ti::String, filename::String)
    melted = melt(top_runs, [:ea, :chromosome])
    top_params = @from i in melted begin
        @where ((i.ea == "oneplus" && i.variable in ea_params) ||
                (i.ea == "GA" && i.variable in ga_params) ||
                (i.ea == "NEAT" && i.variable in neat_params))
        @select i
        @collect DataFrame
    end
    top_params[:parameter] = indexin(top_params[:variable], all_params)

    plt = plot(top_params, xgroup=:ea, ygroup=:chromosome,
               x=:parameter, y=:value, color=:variable,
               Geom.subplot_grid(Geom.boxplot,
               Coord.cartesian(ymin=0.0, ymax=1.0, xmin=0, xmax=length(all_params))),
               Guide.ylabel(nothing), Guide.xlabel(nothing),
               style(key_position=:right),
               Guide.title(ti))
    draw(PDF(string(filename, "_params.pdf"), 12inch, 8inch), plt)
    top_params
end

function get_bests(logs::Array{String})
    bests = DataFrame()
    for l in logs
        res, top_runs = get_cmaes_results(string("results/cmaes/", l, ".log"))
        tops = @from i in top_runs begin
            @where i.id <= 5
            @select i
            @collect DataFrame
        end
        tops[:problem] = l
        bests = vcat(bests, tops)
    end
    bests
end
