using Base.Test
using CGP

const global D = 5

# files = ["cfg/test.yaml", "cfg/classic.yaml", "cfg/matrix.yaml", "cfg/mtcgp.yaml"]
files = ["cfg/test.yaml", "cfg/mtcgp.yaml"]

@testset "Functions" begin
    for file in files
        CGP.Config.init(file)
        @testset "Bounds" begin
            for f in CGP.Config.functions
                println(string(file, ": ", f))
                out = f(rand(), rand(), rand())
                @test all(out .>= -1.0) && all(out .<= 1.0)
                out = f(rand(), rand(rand(1:D, rand(1:D))...), rand())
                @test all(out .>= -1.0) && all(out .<= 1.0)
                out = f(rand(rand(1:D, rand(1:D))...), rand(), rand())
                @test all(out .>= -1.0) && all(out .<= 1.0)
                out = f(rand(rand(1:D, rand(1:D))...), rand(rand(1:D, rand(1:D))...), rand())
                @test all(out .>= -1.0) && all(out .<= 1.0)
            end
        end
        CGP.Config.reset()
    end
end
