using Base.Test
using CGP

# @testset "test functions" begin
#     CGP.Config.init("cfg/test.yaml")
#     for f in CGP.Config.functions
#         out = f(rand(), rand(), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#         out = f(rand(), rand(rand(1:5, rand(1:5))...), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#         out = f(rand(rand(1:5, rand(1:5))...), rand(), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#         out = f(rand(rand(1:5, rand(1:5))...), rand(rand(1:5, rand(1:5))...), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#     end
# end

# @testset "classic functions" begin
#     CGP.Config.init("cfg/classic.yaml")
#     for f in CGP.Config.functions
#         out = f(rand(), rand(), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#         out = f(rand(), rand(rand(1:5, rand(1:5))...), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#         out = f(rand(rand(1:5, rand(1:5))...), rand(), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#         out = f(rand(rand(1:5, rand(1:5))...), rand(rand(1:5, rand(1:5))...), rand())
#         @test all(out .>= 0.0) && all(out .<= 1.0)
#     end
# end

@testset "matrix functions" begin
    CGP.Config.init("cfg/matrix.yaml")
    for f in CGP.Config.functions
        println(repr(f))
        out = f(rand(), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(), rand(rand(1:5), rand(1:5)), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:5), rand(1:5)), rand(), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
        out = f(rand(rand(1:5), rand(1:5)), rand(rand(1:5), rand(1:5)), rand())
        @test all(out .>= 0.0) && all(out .<= 1.0)
    end
end
