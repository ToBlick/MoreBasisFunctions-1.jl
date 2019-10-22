
using BasisFunctions: hasderivative, hasantiderivative

P = 80

@testset "$(rpad("Bernstein polynomials",P))" begin

    ξ = [0.0, 1.0]
    b = Bernstein(ξ)

    @testset "$(rpad("Basic functionality",P))" begin
        @test nodes(b) == ξ
        @test nnodes(b) == 2
        @test degree(b) == 1
        @test hasderivative(b) == true
        @test hasantiderivative(b) == false
        @test support(b) == BernsteinInterval{eltype(ξ)}()
    end


    @testset "$(rpad("Basis evaluation",P))" begin
        @test b[1](0.0) == 1.0
        @test b[1](0.5) == 0.5
        @test b[1](1.0) == 0.0

        @test b[2](0.0) == 0.0
        @test b[2](0.5) == 0.5
        @test b[2](1.0) == 1.0
    end


    @testset "$(rpad("Basis derivative",P))" begin
        @test eval_element_derivative(b, 1, 0.0) == -1.0
        @test eval_element_derivative(b, 1, 0.5) == -1.0
        @test eval_element_derivative(b, 1, 1.0) == -1.0

        @test eval_element_derivative(b, 2, 0.0) == +1.0
        @test eval_element_derivative(b, 2, 0.5) == +1.0
        @test eval_element_derivative(b, 2, 1.0) == +1.0
    end

end
