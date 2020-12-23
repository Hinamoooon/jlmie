using jlmie
using Test

@testset "jlmie" begin
    @test jlmie_Qsca(1,1) == [0.0]
    @test abs.(jlmie_Qsca(4,1) - [6.062]) < [0.01]
end