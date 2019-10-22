module MoreBasisFunctions

    using GridArrays: ScatteredGrid
    using BasisFunctions

    export PolynomIndex
    export Lagrange, LagrangeInterval

    export eval_element, eval_element_derivative, eval_element_antiderivative
    export nodes, nnodes, degree, support, moment

    include("bases/generic/dictionary.jl")
    include("bases/generic/vandermonde_matrix.jl")

    include("bases/poly/polynomials.jl")
    include("bases/poly/lagrange.jl")

end # module
