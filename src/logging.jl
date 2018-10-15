export log_best, log_all

function log_best(id::String, seed::Int64, eval_count::Int64, max_fit::Float64,
                  best::Chromosome, ea::Function, ctype::DataType, log_gen::Bool)
    if log_gen
        logstr = @sprintf("E: %s %d %d %0.5f %d %d",
                          id, seed, eval_count, max_fit,
                          sum([n.active for n in best.nodes]),
                          length(best.nodes))
        if Config.log_config
            logstr = string(logstr,
                            @sprintf("%s %s %s",
                                     string(ea), string(ctype), Config.to_string()))
        end
        @info(logstr)
    end
    if Config.save_best
        @info(@sprintf("C: %s", string(best.genes)))
    end
    nothing
end

function log_all(id::String, seed::Int64, eval_count::Int64, max_fit::Float64,
                  best::Chromosome, ea::Function, ctype::DataType, log_gen::Bool)
    log_best(id, seed, eval_count, max_fit, best, ea, ctype, true)
end
