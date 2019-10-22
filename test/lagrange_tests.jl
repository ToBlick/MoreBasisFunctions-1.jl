
using MoreBasisFunctions: eval_element_antiderivative
using BasisFunctions: eval_element_derivative, moment
using BasisFunctions: hasderivative, hasantiderivative, support
using BasisFunctions: PolynomialDegree

ξ = [-1.0, +1.0]
b = Lagrange(ξ)
P = 80

@testset "$(rpad("Basic functionality",P))" begin
    @test nodes(b) == ξ
    @test nnodes(b) == 2
    @test degree(b) == 1
    @test hasderivative(b) == true
    @test hasantiderivative(b) == true
    @test support(b) == LagrangeInterval{eltype(ξ)}()
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
    @test moment(b, 1) == 1.0
    @test moment(b, 2) == 1.0

    @test eval_element_antiderivative(b, 1, -1.0) == 0.0
    @test eval_element_antiderivative(b, 1,  0.0) == 0.75
    @test eval_element_antiderivative(b, 1, +1.0) == 1.0

    @test eval_element_antiderivative(b, 2, -1.0) == 0.0
    @test eval_element_antiderivative(b, 2,  0.0) == 0.25
    @test eval_element_antiderivative(b, 2, +1.0) == 1.0
end
