module MoreBasisFunctions

    using GridArrays: ScatteredGrid
    using BasisFunctions
    using BasisFunctions: PolynomialBasis, PolynomialDegree

    export PolynomIndex
    export Lagrange, LagrangeInterval

    export nodes, nnodes, degree

    include("bases/generic/dictionary.jl")
    include("bases/generic/vandermonde_matrix.jl")

    const PolynomIndex = Int

    include("bases/poly/lagrange.jl")


end # module
