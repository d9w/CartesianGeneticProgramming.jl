using Gadfly
using Distributions
using Colors
using DataFrames
using Query
using ColorSchemes

Gadfly.push_theme(Theme(major_label_font="Helvetica",
                        minor_label_font="Helvetica",
                        major_label_font_size=18pt, minor_label_font_size=16pt,
                        line_width=0.8mm, key_label_font="Droid Sans",
                        point_size=1.0mm,
                        lowlight_color=c->RGBA{Float32}(c.r, c.g, c.b, 0.2),
                        key_label_font_size=14pt,
                        default_color=colorant"#000000"))

colors = [colorant"#e41a1c", colorant"#377eb8", colorant"#4daf4a",
          colorant"#984ea3", colorant"#ff7f00", colorant"#ffff33"]

pclasses = ["classification", "regression", "rl"]
shortclass = ["Classification", "Regression", "Reinforcement Learning"]

mutations = ["gene_mutate", "mixed_node_mutate", "mixed_subtree_mutate"]

crossovers = ["single_point_crossover", "proportional_crossover",
              "random_node_crossover", "aligned_node_crossover",
              "output_graph_crossover", "subgraph_crossover"]

function stdz(x)
    y = std(x)
    if isnan(y)
        y = 0.0
    end
    y
end

function get_param_results(log::String)
    res = readtable(log, header=false, separator=' ',
                    names=[:id, :problem, :pclass, :seed, :eval, :fit, :refit,
                           :nactive, :nodes, :ea, :chromo, :mut, :active,
                           :cross, :weights, :node_min, :node_static, :node_max,
                           :input_start, :recurrency, :input_mutation_rate,
                           :output_mutation_rate, :node_mutation_rate,
                           :node_size_delta, :modify_mutation_rate, :lambda,
                           :ga_population, :ga_elitism_rate, :ga_crossover_rate,
                           :ga_mutation_rate])
    res[:ea][res[:ea].=="oneplus"] = "1+λ"
    res[:chromosome] = map(x->String(split(split(x, '.')[2], "Chromo")[1]),
                                   res[:chromo])
    res[:mutation] = map(x->(findfirst(x.==mutations)-1)/(length(mutations)-1),
                                 res[:mut])
    res[:crossover] = map(x->(findfirst(x.==crossovers)-1)/(length(crossovers)-1),
                                  res[:cross])
    res[:class] = map(x->shortclass[(findfirst(x.==pclasses))], res[:pclass])
    res[:active] = Int64.(res[:active])
    res[:weights] = Int64.(res[:weights])
    res[:lambda] = res[:lambda] / 10
    res[:ga_population] = res[:ga_population] / 200
    res[:input_start] = 1.0.+res[:input_start]
    gene_mutate_ids = res[:mut] .== "gene_mutate"
    res[:modify_mutation_rate][gene_mutate_ids] = NaN
    res[:node_size_delta][gene_mutate_ids] = NaN
    res[:uuid] = string.(res[:id], "_", res[:pclass], "_", res[:ea], "_",
                         res[:chromosome])
    res
end

function get_top_by_pareto(res::DataFrame, nbests::Int64=10)
    full_ranks = DataFrame()
    for subdf in groupby(res, :pclass)
        pmeans = by(subdf, [:uuid, :problem], df->DataFrame(mfit=mean(df[:fit])))
        if split(pmeans[:uuid][1], "_")[2] == "classification"
            pmeans[:prob_id] = string.("p", map(x->findfirst(
                x.==["cancer", "diabetes", "glass"]), pmeans[:problem]))
        elseif split(pmeans[:uuid][1], "_")[2] == "regression"
            pmeans[:prob_id] = string.("p", map(x->findfirst(
                x.==["abalone", "fires", "quality"]), pmeans[:problem]))
        else
            pmeans[:prob_id] = string.("p", map(x->findfirst(
                x.==["ant", "cheetah", "humanoid"]), pmeans[:problem]))
        end
        fits = unstack(pmeans, :uuid, :prob_id, :mfit)
        fits = fits[((!isna.(fits[:p1])).&(!isna.(fits[:p2])).&(!isna.(fits[:p3]))), :]
        fits[:ea] = map(x->split(x, "_")[3], fits[:uuid])
        fits[:chromosome] = map(x->split(x, "_")[4], fits[:uuid])
        for subfits_ea in groupby(fits, :ea)
            for scres in groupby(subfits_ea, :chromosome)
                ranks = DataFrame()
                ranks[:uuid] = scres[:uuid]
                # ranks[:score] = map(
                #     x->sum((scres[:p1][x].>scres[:p1]).&(scres[:p2][x].>scres[:p2]).&
                #           (scres[:p3][x].>scres[:p3])), 1:size(scres, 1))
                ranks[:score] = (scres[:p1] + scres[:p2] + scres[:p3]) / 3.0
                sort!(ranks, cols=(order(:score, rev=true)))
                rranks = head(ranks, nbests)
                maxscore = maximum(rranks[:score])
                minscore = minimum(rranks[:score])
                rranks[:rank] = 1:size(rranks[:score], 1)
                rranks[:nscore] = (rranks[:score] - minscore) / (maxscore - minscore)
                full_ranks = vcat(full_ranks, rranks)
            end
        end
    end
    joined = join(full_ranks, res, on=:uuid, kind=:inner)
    best_params = aggregate(joined, :uuid, first)
    names!(best_params, names(joined))
    sort!(best_params, cols=(:pclass, :ea, :chromosome, :rank))
    best_params
end

function get_norm_means(res::DataFrame)
    problem_stats = by(res, [:problem], df->DataFrame(
        pmax=maximum(df[:fit]), pmin=minimum(df[:fit])))

    norm_fit = zeros(size(res,1))
    for problem in unique(problem_stats[:problem])
        stat_id = findfirst(problem_stats[:problem].==problem)
        pmax = problem_stats[:pmax][stat_id]
        pmin = problem_stats[:pmin][stat_id]
        println(problem, " ", pmin, " ", pmax)
        pinds = res[:problem] .== problem
        norm_fit[pinds] = (res[:fit][pinds] - pmin)./(pmax - pmin)
    end

    res[:norm_fit] = norm_fit

    full_nmeans = by(res, [:id, :pclass, :ea, :chromosome], df->DataFrame(
        mfit=mean(df[:norm_fit]), sfit=std(df[:norm_fit]),
        nprob=length(unique(df[:problem])), mut=df[:mut][1],
        active=df[:active][1], cross=df[:cross][1], weights=df[:weights][1],
        input_start=df[:input_start][1], recurrency=df[:recurrency][1],
        input_mutation_rate=df[:input_mutation_rate][1],
        output_mutation_rate=df[:output_mutation_rate][1],
        node_mutation_rate=df[:node_mutation_rate][1],
        node_size_delta=df[:node_size_delta][1],
        modify_mutation_rate=df[:modify_mutation_rate][1],
        lambda=df[:lambda][1], ga_population=df[:ga_population][1],
        ga_elitism_rate=df[:ga_elitism_rate][1],
        ga_crossover_rate=df[:ga_crossover_rate][1],
        ga_mutation_rate=df[:ga_mutation_rate][1]))

    nmeans = @from i in full_nmeans begin
        @where i.nprob == 3
        @select i
        @collect DataFrame
    end

    nmeans
end

function get_best_params(nmeans::DataFrame, top::Int64=20)
    best_params = DataFrame()
    # best 20 per pclass, ea, chromo
    # 20 * 3 * 2 * 2
    for class in unique(nmeans[:pclass])
        for ea in unique(nmeans[:ea])
            for chromosome in unique(nmeans[:chromosome])
                class_res = @from i in nmeans begin
                    @where ((i.pclass == class) && (i.ea==ea) &&
                            (i.chromosome == chromosome))
                    @select i
                    @collect DataFrame
                end
                sort!(class_res, cols=(order(:mfit, rev=true)))
                best_params = vcat(best_params, head(class_res, top))
            end
        end
    end
    best_params
end

function plot_fits(best_params::DataFrame)
    plts = []
    for class in shortclass
        class_params = @from i in best_params begin
            @where i.class == class
            @select i
            @collect DataFrame
        end
        plt = plot(class_params, xgroup=:ea, x=:chromosome, y=:mfit,
                   Guide.xlabel(nothing), Guide.ylabel(nothing),
                   Guide.title(class),
                   Geom.subplot_grid(Geom.violin))
        append!(plts, [plt])
    end
    final_plt = hstack(plts...)
    draw(PDF("best_fits.pdf", 12inch, 4inch), final_plt)
end

function rename_melted!(melted::DataFrame)
    names = Array{String}(size(melted, 1))
    names[melted[:variable].==:mfit] = "fitness"
    names[melted[:variable].==:mutation] = "mutation"
    names[melted[:variable].==:crossover] = "crossover"
    names[melted[:variable].==:lambda] = "population"
    names[melted[:variable].==:ga_population] = "population"
    names[melted[:variable].==:weights] = "w"
    names[melted[:variable].==:recurrency] = "r"
    names[melted[:variable].==:active] = "m<sub>active</sub>"
    names[melted[:variable].==:input_start] = "I<sub>start</sub>"
    names[melted[:variable].==:input_mutation_rate] = "m<sub>input</sub>"
    names[melted[:variable].==:output_mutation_rate] = "m<sub>output</sub>"
    names[melted[:variable].==:node_mutation_rate] = "m<sub>node</sub>"
    names[melted[:variable].==:node_size_delta] = "m<sub>δ</sub>"
    names[melted[:variable].==:modify_mutation_rate] = "m<sub>modify</sub>"
    names[melted[:variable].==:ga_elitism_rate] = "GA<sub>elitism</sub>"
    names[melted[:variable].==:ga_crossover_rate] = "GA<sub>crossover</sub>"
    names[melted[:variable].==:ga_mutation_rate] = "GA<sub>mutation</sub>"
    melted[:names] = names
    melted
end

function get_filename(ea::String, chromo::String)
    sea = ea
    if ea == "1+λ"
        sea = "oneplus"
    end
    string(sea, "_", chromo, "_params.pdf")
end

function get_title(pclass::String)
    sclass = string(uppercase(pclass[1]), pclass[2:end])
    if pclass == "rl"
        sclass = "Reinforcement Learning"
    end
    sclass
end

function to_yaml(best_params::DataFrame, class::String, ea::String,
                 chromosome::String)
    ps = @from i in best_params begin;
        @where ((i.class==class)&(i.ea==ea)&(i.chromosome==chromosome))
        @select i; @collect DataFrame
    end
    melted = melt(ps, [:uuid])
    cols = [:mut, :cross, :lambda, :ga_population, :input_start,
            :node_min, :node_static, :node_max, :eval, :recurrency, :weights,
            :active, :input_mutation_rate, :output_mutation_rate,
            :node_mutation_rate, :node_size_delta, :modify_mutation_rate,
            :ga_elitism_rate, :ga_crossover_rate, :ga_mutation_rate]
    cfg = ""
    # cfg = Dict([("save_best", false), ("node_inputs", 2)])
    cfg = string(cfg, "save_best: false\nnode_inputs: 2")
    for i in 1:size(melted, 1)
        var = melted[:variable][i]
        val = melted[:value][i]
        if var in cols
            if var == :node_min
                cfg = string(cfg, "\nstarting_nodes: ", val)
            elseif var == :node_static
                cfg = string(cfg, "\nstatic_node_size: ", val)
            elseif var == :node_max
                cfg = string(cfg, "\nnode_size_cap: ", val)
            elseif var == :eval
                cfg = string(cfg, "\ntotal_evals: ", val)
            elseif var == :mut
                cfg = string(cfg, "\nmutate_method: \":", val, "\"")
            elseif var == :cross
                cfg = string(cfg, "\ncrossover_method: \":", val, "\"")
            elseif var == :active
                cfg = string(cfg, "\nactive_mutate: ", Bool(val))
            elseif var == :weights
                cfg = string(cfg, "\nweights: ", Bool(val))
            elseif var == :lambda
                cfg = string(cfg, "\nlambda: ", Int64(10*val))
            elseif var == :ga_population
                cfg = string(cfg, "\nga_population: ", Int64(200*val))
            elseif var == :input_start
                cfg = string(cfg, "\ninput_start: ", val-1.0)
            else
                if isnan(val); val = 0.0; end
                cfg = string(cfg, "\n", var, ": ", val)
            end
        end
    end
    string(cfg, "\n")
end

function get_correlation(res::DataFrame)
    all_cors = DataFrame()
    cols = [:mutation, :crossover, :lambda, :ga_population, :input_start,
            :recurrency, :weights, :active, :input_mutation_rate,
            :output_mutation_rate, :node_mutation_rate, :node_size_delta,
            :modify_mutation_rate, :ga_elitism_rate, :ga_crossover_rate,
            :ga_mutation_rate]
    pmeans = by(res, [:uuid, :problem], df->DataFrame(mfit=mean(df[:fit])))
    joined = join(pmeans, res, on=:uuid, kind=:inner)
    for c in cols
        joined[c][isnan(joined[c])] = 0.0
        cdf = by(joined, [:ea, :chromosome, :pclass, :class],
                 df->DataFrame(cors=abs(cor(df[:mfit], df[c])), variable=c))
        cdf[:cors][isnan(cdf[:cors])] = 0.0
        all_cors = vcat(all_cors, cdf)
    end
    all_cors = rename_melted!(all_cors)
    all_cors[:method] = string.(all_cors[:ea], " ", all_cors[:chromosome])
    sort!(all_cors, cols=(:class, :ea, :chromosome))
    all_cors
end

function plot_correlations(cors::DataFrame)
    plts = []
    methods = unique(cors[:method])
    cl = sortcolorscheme(ColorSchemes.fuchsia, rev=true)
    for method in methods
        cres = @from i in cors begin
            @where (i.method==method)
            @select i
            @collect DataFrame
        end
        sort!(cres, cols=(order(:pclass, rev=true)))
        maxcor = maximum(cres[:cors])
        # mincor = minimum(cres[:cors])
        cres[:cor] = maxcor - cres[:cors]
        # println(cres[:cor])
        # cres[:cor][isinf(cres[:cor])] = 100.0
        plt = plot(cres, x=:names, y=:class, color=:cors, Geom.rectbin,
                   Guide.colorkey(title="cor"),
                   Scale.ContinuousColorScale(p -> get(cl, p)),
                   Guide.xlabel(nothing), Guide.ylabel(nothing),
                   Guide.title(method));
        append!(plts, [plt])
    end
    plt = vstack(plts...)
    draw(PDF("correlations.pdf", 20inch, 10inch), plt)
    nothing
end

function plot_params(bests::DataFrame, ea::String, chromosome::String, cols::Array{Symbol})
    plts = []
    cres = @from i in bests begin
        @where ((i.ea==ea) && (i.chromosome==chromosome))
        @select i
        @collect DataFrame
    end
    # cres[:rank] = Int64.(repmat(1:(size(cres,1)/3), 3))
    melted = melt(cres, [:nscore, :class], cols)
    melted = rename_melted!(melted)
    plt = plot(melted, x=:names, ygroup=:class, y=:value, color=:nscore,
               Geom.subplot_grid(Geom.beeswarm),
               Guide.xlabel(nothing), Guide.ylabel(nothing),
               Guide.colorkey(title="score"),
                # Guide.xticks(orientation=:vertical),
               Scale.ContinuousColorScale(p -> get(ColorSchemes.fuchsia, 1.0-p)),
               # Scale.color_discrete_manual(colors...),
               Guide.title(string(ea, " ", chromosome)))
    draw(PDF(get_filename(ea, chromosome), 24inch, 10inch), plt)
    nothing
end

function plot_all(bests::DataFrame)
    plot_params(bests, "1+λ", "CGP",
                [:mutation, :lambda, :recurrency, :weights,
                 :active, :output_mutation_rate, :node_mutation_rate,
                 :node_size_delta, :modify_mutation_rate]);
    plot_params(bests, "1+λ", "PCGP",
                [:mutation, :lambda, :input_start, :recurrency, :weights,
                 :active, :input_mutation_rate, :output_mutation_rate, :node_mutation_rate,
                 :node_size_delta, :modify_mutation_rate])
    plot_params(bests, "GA", "CGP",
                [:mutation, :crossover, :ga_population, :recurrency, :weights,
                 :active, :output_mutation_rate, :node_mutation_rate,
                 :node_size_delta, :modify_mutation_rate,
                 :ga_elitism_rate, :ga_crossover_rate, :ga_mutation_rate])
    plot_params(bests, "GA", "PCGP",
                [:mutation, :crossover, :ga_population, :input_start, :recurrency,
                 :weights, :active, :input_mutation_rate, :output_mutation_rate,
                 :node_mutation_rate, :node_size_delta, :modify_mutation_rate,
                 :ga_elitism_rate, :ga_crossover_rate, :ga_mutation_rate])
end
