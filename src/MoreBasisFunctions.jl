module MoreBasisFunctions

    using GridArrays: ScatteredGrid
    using BasisFunctions
    using BasisFunctions: PolynomialBasis, PolynomialDegree

    export Lagrange, LagrangeInterval

    export nodes, nnodes, degree
    export hasderivative, hasantiderivative

    include("bases/general/vandermonde_matrix.jl")

    unsafe_eval_element_antiderivative(dict::Dictionary, idx, x) =
        unsafe_eval_element_antiderivative(dict, native_index(dict, idx), x)

    include("bases/poly/lagrange.jl")


end # module
