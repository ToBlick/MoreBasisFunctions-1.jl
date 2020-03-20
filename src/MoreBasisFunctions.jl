module MoreBasisFunctions

    using GridArrays: ScatteredGrid
    using BasisFunctions

    export PBSpline, PBSplineIndex, PBSplineInterval
    export Bernstein, BernsteinIndex, BernsteinInterval
    export Lagrange, LagrangeIndex, LagrangeInterval
    export LagrangeGaussLegendre, LagrangeLobattoLegendre
    export LagrangeGaussChebyshev, LagrangeLobattoChebyshev

    export eval_element, eval_element_derivative, eval_element_antiderivative
    export nodes, nnodes, degree, support

    include("bases/generic/dictionary.jl")
    include("bases/generic/nodes.jl")
    include("bases/generic/vandermonde_matrix.jl")

    include("bases/poly/bernstein.jl")
    include("bases/poly/lagrange.jl")
    include("bases/poly/periodicbspline.jl")

end # module
