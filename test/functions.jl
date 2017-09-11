using Base.Test
using CGP

const global D = 5

# files = ["cfg/test.yaml", "cfg/classic.yaml", "cfg/matrix.yaml", "cfg/mtcgp.yaml"]
files = ["cfg/test.yaml", "cfg/mtcgp.yaml"]

function testout(out)
    minout = minimum(out); maxout = maximum(out)
    @printf("%0.2f, %0.2f, ", minout, maxout)
    @test (minout >= -1.0) && (maxout <= 1.0)
end


@testset "Functions" begin
    for file in files
        CGP.Config.init(file)
        @testset "$file bounds" begin
            for f in CGP.Config.functions
                print(string(file, ": ", f, ", "))
                out = f(rand(), rand(), rand())
                testout(out)
                out = f(rand(), rand(rand(1:D, rand(1:D))...), rand())
                testout(out)
                out = f(rand(rand(1:D, rand(1:D))...), rand(), rand())
                testout(out)
                out = f(rand(rand(1:D, rand(1:D))...), rand(rand(1:D, rand(1:D))...), rand())
                testout(out)
                println()
            end
        end
        CGP.Config.reset()
    end
end
