using BasisFunctions: hasderivative, hasantiderivative, resize

P = 80

b = PBSpline(2, 5)

@testset "$(rpad("Periodic BSplines",P))" begin

    @testset "$(rpad("Basic functionality",P))" begin
        ξ = nodes(b)
        @test nnodes(b) == 5
        @test degree(b) == 2
        @test support(b) == PBSplineInterval{eltype(ξ)}()
        @test nodes(b) == nodes(similar(b, eltype(ξ), 5))
        @test nodes(b) == nodes(similar(b, eltype(ξ), ξ))
        @test nodes(resize(b, 5)) == ξ
        @test nodes(resize(b, 3)) == [0.0, 0.5, 1.0]
    end

    @testset "$(rpad("Basis evaluation",P))" begin
        @test b[1](0.2) == 0.5
        @test b[1](0.3) == 0.75
        @test b[1](0.6) == 0.0

        @test b[5](0.0) == 0.5
        @test b[5](1.0) == 0.5
        @test b[5](0.9) == 0.125
    end

    @testset "$(rpad("Basis derivative",P))" begin
        @test eval_element_derivative(b, 1, 0.1) == 2.5
        @test eval_element_derivative(b, 1, 0.3) == 0.0
    end

end
