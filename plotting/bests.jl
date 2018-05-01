using Gadfly
using Distributions
using Colors
using DataFrames
using Query
using HypothesisTests
using LaTeXStrings

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

function reducef(df, xmax)
    r = 1:length(df[:eval])
    if df[:eval][end] != xmax
        r = 1:(length(df[:eval])+1)
    end
    map(i->mapf(i, df, xmax), r)
end

function mapf(i::Int64, df, xmax::Int64)
    if i > length(df[:eval])
        if df[:eval][end] < xmax
            return df[:fit][end] * ones(xmax - df[:eval][end])
        else
            return []
        end
    end
    lower = 0
    if i >= 2
        lower = df[:eval][i-1]
    end
    if df[:eval][i] >= lower
        return df[:fit][i] * ones(df[:eval][i] - lower)
    else
        println(df[:ea][i], " ", df[:chromosome][i], " ", df[:eval][i], " ",
                df[:seed][i], " ", lower)
        return []
    end
end

function get_problem_stats(res::DataFrame, problem::String;
                           xmax::Int64=maximum(res[:eval]))
    pres = @from i in res begin
        @where i.problem == problem
        @select i; @collect DataFrame;
    end
    filled = by(pres, [:seed, :id],
                df->reduce(vcat, reducef(df, xmax)))
    filled[:xs] = repeat(1:xmax, outer=Int64(size(filled,1)/xmax))
    stats = by(filled, [:xs, :id],
               df->DataFrame(stds=std(df[:x1]), means=mean(df[:x1]),
                             mins=minimum(df[:x1]),
                             maxs=maximum(df[:x1])))
    stats[:stds][isnan.(stats[:stds])] = 0;
    stats[:lower] = stats[:means]-0.5*stats[:stds]
    stats[:upper] = stats[:means]+0.5*stats[:stds]
    stats
end

function to_id(ea::String, chromosome::String)
    id = "e<sub>1</sub>"
    if ea == "1+λ"
        if chromosome == "PCGP"
            id = "e<sub>2</sub>"
        end
    else
        if chromosome == "CGP"
            id = "e<sub>3</sub>"
        else
            id = "e<sub>4</sub>"
        end
    end
    id
end

function get_bests_results(log::String)
    res = readtable(log, header=false, separator=' ',
                    names=[:problem, :pclass, :seed, :eval, :fit, :refit,
                           :nactive, :nodes, :ea, :chromo])
    res[:ea][res[:ea].=="oneplus"] = "1+λ"
    res[:chromosome] = map(x->String(split(split(x, '.')[2], "Chromo")[1]),
                                   res[:chromo])
    res[:id] = to_id.(res[:ea], res[:chromosome])
    res[:method] = string.(res[:ea], " ", res[:chromosome])
    sort!(res, cols=[:id])
    res
end

function plot_evolution(stats::DataFrame, title::String,
                        filename::String="training";
                        xmin=minimum(stats[:xs]),
                        xmax=maximum(stats[:xs]),
                        ymin=minimum(stats[:lower]),
                        ymax=maximum(stats[:upper]),
                        ylabel="Fitness",
                        key_position=:right)
    plt = plot(stats, x="xs", y="means", ymin="lower", ymax="upper", color="id",
               Geom.line, Geom.ribbon,
               Scale.color_discrete_manual(colors...),
               Guide.title(title),
               Guide.xlabel("Evaluations"),
               Guide.ylabel(ylabel),
               Coord.cartesian(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
               style(key_position=key_position))
    draw(PDF(string("plots/", filename, ".pdf"), 8inch, 6inch), plt)
    plt
end

function plot_bests(res::DataFrame, pclass::String, tstr::String,
                    problems::Array{String}, labels::Array{String};
                    xmax::Int64=maximum(res[:eval]))
    plts = []
    key_positions = [:none, :none, :right]
    for p in eachindex(problems)
        cres = @from i in res begin
            @where ((i.problem==problems[p])&(i.eval==xmax))
            @select i; @collect DataFrame;
        end

        plt = plot(cres, x=:id, y=:refit, color=:id,
                   Geom.boxplot,
                   Scale.color_discrete_manual(colors...),
                   Guide.title(labels[p]),
                   Guide.xlabel(nothing),
                   Guide.ylabel(nothing),
                   style(key_position=key_positions[p]))
        append!(plts, [plt])
    end

    plt = hstack(plts...)
    plt = title(plt, tstr)
    draw(PDF(string("plots/", pclass, "_bests.pdf"), 18inch, 4inch), plt)
    nothing
end

function pvalues(res::DataFrame, pclass::String; xmax::Int64=maximum(res[:eval]))
    cres = @from i in res begin
        @where ((i.pclass==pclass)&(i.eval==xmax))
        @select i; @collect DataFrame;
    end
    pv = DataFrame()
    methods = unique(res[:method])
    problems = Array{String}(unique(cres[:problem]))
    for p in eachindex(problems)
        pres = @from i in cres begin
            @where (i.problem==problems[p])
            @select i; @collect DataFrame;
        end
        s = string(problems[p], "\n& ")
        for n in eachindex(methods)
            s = string(s, L"$e_", n, L"$", " & ")
        end
        s = string(s, "\n\\hline\n")
        for n in eachindex(methods)
            s = string(s, L"$e_", n, L"$", " & ")
            nres = @from i in pres begin
                @where (i.method==methods[n])
                @select i.refit; @collect
            end
            nres = Array{Float64}(nres)
            for m in eachindex(methods)
                mres = @from i in pres begin
                    @where (i.method==methods[m])
                    @select i.refit; @collect
                end
                mres = Array{Float64}(mres)
                # println(methods[n], methods[m], pvalue(OneSampleTTest(mres, nres)))
                pval = pvalue(OneSampleTTest(mres, nres))
                sval = @sprintf("%0.3f", pval)
                if isnan(pval)
                    pval = 0.0
                    sval = " "
                end
                s = string(s, sval, " & ")
                pv = vcat(pv, DataFrame(problem=problems[p], m1=methods[n],
                                        m2=methods[m], p=pval))
            end
            s = string(s, "\n")
        end
        s = replace(s, "1+λ", L"$1+\lambda$")
        s = replace(s, " & \n", "\\\\\n")
        println(s)
    end
    pv
end
