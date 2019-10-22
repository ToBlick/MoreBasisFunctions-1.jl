
using BasisFunctions: hasderivative, hasantiderivative

P = 80

@testset "$(rpad("Lagrange polynomials",P))" begin

    ξ = [-1.0, +1.0]
    b = Lagrange(ξ)
    blob = LagrangeLobatto(2)
    bleg = LagrangeLegendre(2)


    @testset "$(rpad("Basic functionality",P))" begin
        @test nodes(b) == ξ
        @test nnodes(b) == 2
        @test degree(b) == 1
        @test hasderivative(b) == true
        @test hasantiderivative(b) == true
        @test support(b) == LagrangeInterval{eltype(ξ)}()
        @test nodes(blob) == nodes(b)
        @test nodes(bleg) == [-0.5773502691896258, 0.5773502691896258]
        @test nodes(blob) == nodes(similar(blob, eltype(ξ), 2))
        @test nodes(bleg) == nodes(similar(bleg, eltype(ξ), 2))
    end


    @testset "$(rpad("Basis evaluation",P))" begin
        @test b[1](-1.0) == 1.0
        @test b[1]( 0.0) == 0.5
        @test b[1](+1.0) == 0.0

        @test b[2](-1.0) == 0.0
        @test b[2]( 0.0) == 0.5
        @test b[2](+1.0) == 1.0
    end


    @testset "$(rpad("Basis derivative",P))" begin
        @test eval_element_derivative(b, 1, -1.0) == -0.5
        @test eval_element_derivative(b, 1,  0.0) == -0.5
        @test eval_element_derivative(b, 1, +1.0) == -0.5

        @test eval_element_derivative(b, 2, -1.0) == +0.5
        @test eval_element_derivative(b, 2,  0.0) == +0.5
        @test eval_element_derivative(b, 2, +1.0) == +0.5
    end


    @testset "$(rpad("Basis antiderivative",P))" begin

        @test eval_element_antiderivative(b, 1, -1.0) == 0.0
        @test eval_element_antiderivative(b, 1,  0.0) == 0.75
        @test eval_element_antiderivative(b, 1, +1.0) == 1.0

        @test eval_element_antiderivative(b, 2, -1.0) == 0.0
        @test eval_element_antiderivative(b, 2,  0.0) == 0.25
        @test eval_element_antiderivative(b, 2, +1.0) == 1.0
    end

end
