using Base.Test
using CGP

const global D = 5

files = ["cfg/test.yaml", "cfg/classic.yaml", "cfg/mtcgp.yaml"]

function test_inps(f, inps)
    out = f(inps...)
    @test all(out == f(inps...))
    minout = minimum(out); maxout = maximum(out)
    @printf("%0.2f, %0.2f, %s, ", minout, maxout, typeof(out))
    @test (minout >= -1.0) && (maxout <= 1.0)
end


@testset "Functions" begin
    for file in files
        CGP.Config.init(file)
        @testset "$file" begin
            for f in CGP.Config.functions
                print(string(file, ": ", f, ", "))
                test_inps(f, [rand(), rand(), rand()])
                test_inps(f, [rand(), rand(rand(1:D, rand(1:D))...), rand()])
                test_inps(f, [rand(rand(1:D, rand(1:D))...), rand(), rand()])
                test_inps(f, [rand(rand(1:D, rand(1:D))...), rand(rand(1:D, rand(1:D))...),
                           rand()])
                println()
            end
        end
        CGP.Config.reset()
    end
end
