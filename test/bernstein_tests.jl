
using BasisFunctions: hasderivative, hasantiderivative, resize

P = 80

@testset "$(rpad("Bernstein polynomials",P))" begin

    ξ = [0.0, 1.0]
    b = Bernstein(ξ)

    BasisFunctions.Test.test_generic_dict_interface(Bernstein(10))

    @testset "$(rpad("Basic functionality",P))" begin
        @test nodes(b) == ξ
        @test nnodes(b) == 2
        @test degree(b) == 1
        @test support(b) == BernsteinInterval{eltype(ξ)}()
        @test nodes(b) == nodes(similar(b, eltype(ξ), 2))
        @test nodes(b) == nodes(similar(b, eltype(ξ), ξ))
        @test nodes(resize(b, 2)) == ξ
        @test nodes(resize(b, 5)) == [0.0, 0.25, 0.50, 0.75, 1.0]
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
