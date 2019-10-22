module MoreBasisFunctions

    using GridArrays: ScatteredGrid
    using BasisFunctions

    export Bernstein, BernsteinIndex, BernsteinInterval
    export Lagrange, LagrangeIndex, LagrangeInterval
    export LagrangeLegendre, LagrangeLobatto

    export eval_element, eval_element_derivative, eval_element_antiderivative
    export nodes, nnodes, degree, support

    include("bases/generic/dictionary.jl")
    include("bases/generic/vandermonde_matrix.jl")

    include("bases/poly/bernstein.jl")
    include("bases/poly/lagrange.jl")

end # module
