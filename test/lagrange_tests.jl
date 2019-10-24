
using BasisFunctions: hasderivative, hasantiderivative

P = 80

@testset "$(rpad("Lagrange polynomials",P))" begin

    ξ = [-1.0, +1.0]
    b = Lagrange(ξ)
    blob = LagrangeLobattoLegendre(2)
    bgau = LagrangeGaussLegendre(2)
    blch = LagrangeLobattoChebyshev(2)
    bgch = LagrangeGaussChebyshev(2)

    BasisFunctions.Test.test_generic_dict_interface(LagrangeGaussLegendre(10))

    @testset "$(rpad("Basic functionality",P))" begin
        @test nodes(b) == ξ
        @test nnodes(b) == 2
        @test degree(b) == 1
        @test support(b) == LagrangeInterval{eltype(ξ)}()
        @test nodes(blob) == nodes(b)
        @test nodes(blch) == nodes(b)
        @test nodes(bgau) == [-0.5773502691896258, 0.5773502691896258]
        @test nodes(bgch) == [-0.7071067811865475, 0.7071067811865475]
        @test nodes(blob) == nodes(similar(blob, eltype(ξ), 2))
        @test nodes(bgau) == nodes(similar(bgau, eltype(ξ), 2))
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
        # @test moment(b, 1) == 1.0
        # @test moment(b, 2) == 1.0

        @test eval_element_antiderivative(b, 1, -1.0) == 0.0
        @test eval_element_antiderivative(b, 1,  0.0) == 0.75
        @test eval_element_antiderivative(b, 1, +1.0) == 1.0

        @test eval_element_antiderivative(b, 2, -1.0) == 0.0
        @test eval_element_antiderivative(b, 2,  0.0) == 0.25
        @test eval_element_antiderivative(b, 2, +1.0) == 1.0
    end

end
