using Base.Test
using CGP

const global D = 5

@testset "test functions" begin
    CGP.Config.init("cfg/test.yaml")
    for f in CGP.Config.functions
        out = f(rand(), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(), rand(rand(1:D, rand(1:D))...), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:D, rand(1:D))...), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:D, rand(1:D))...), rand(rand(1:D, rand(1:D))...), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
    end
end

@testset "classic functions" begin
    CGP.Config.init("cfg/classic.yaml")
    for f in CGP.Config.functions
        out = f(rand(), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(), rand(rand(1:D, rand(1:D))...), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:D, rand(1:D))...), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:D, rand(1:D))...), rand(rand(1:D, rand(1:D))...), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
    end
end

@testset "matrix functions" begin
    CGP.Config.init("cfg/matrix.yaml")
    for f in CGP.Config.functions
        debug(repr(f))
        out = f(rand(), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(), rand(rand(1:D), rand(1:D)), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:D), rand(1:D)), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:D), rand(1:D)), rand(rand(1:D), rand(1:D)), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
    end
end
