using BasisFunctions: UnitInterval, PolynomialBasis
using BasisFunctions: hasderivative, hasantiderivative, support

const PBSplineInterval = UnitInterval
const PBSplineIndex = NativeIndex{:pbspline}

# standard grid generation
get_pbspline_nodes(T,n) = collect(LinRange{T}(0.0, 1.0, n))

function check_equispaced(ξ::Vector{T}) where {T}
    Δ = ξ[2] - ξ[1]
    for i in 1:(length(ξ)-1)
        abs(ξ[i+1] - ξ[i] - Δ) > 1/length(ξ) * 1e-2 ? (return false) : nothing
    end
    true
end

# constructors
struct PBSpline{T} <: PolynomialBasis{T,T}
    n :: Int                    # spline degree
    nodes :: ScatteredGrid{T}   # knot points

    function PBSpline(n::Int, nodes::ScatteredGrid{T}) where {T}
        @assert length(nodes) > n
        @assert check_equispaced(nodes.points)

        new{T}(n, nodes)
    end

    function PBSpline(n::Int, ξ::Vector{T}) where {T}
        PBSpline(n, ScatteredGrid(ξ, PBSplineInterval{T}()))
    end
end

# constructors for standard grid
PBSpline(n::Int) = PBSpline(n-1, get_pbspline_nodes(Float64,n))
PBSpline{T}(n::Int) where {T} = PBSpline(n-1, get_pbspline_nodes(T,n))
PBSpline(ξ::Vector{T}) where {T} = PBSpline(length(ξ)-1, ξ)

# constructor for interval [a,b)
PBSpline(n, ξ, a::Number, b::Number) = rescale(PBSpline(ξ), a, b)

# basis properties
nodes(b::PBSpline)  = b.nodes.points
nnodes(b::PBSpline) = length(b.nodes)
degree(b::PBSpline) = b.n
Base.size(b::PBSpline) = (nnodes(b),)

BasisFunctions.native_index(b::PBSpline, idx) = PBSplineIndex(idx)
BasisFunctions.linear_index(b::PBSpline, idx) = PBSplineIndex(idx)
BasisFunctions.support(b::PBSpline{T}) where {T} = PBSplineInterval{T}()

BasisFunctions.interpolation_grid(b::PBSpline{T}) where {T} = ScatteredGrid(get_pbspline_nodes(T,length(b)), PBSplineInterval{T}())

BasisFunctions.similar(::PBSpline, ::Type{T}, n::Int) where {T} = PBSpline{T}(n)
BasisFunctions.similar(::PBSpline, ::Type{T}, ξ::Vector{T}) where {T} = PBSpline(ξ)

# evaluation at point x
function _pbspline(b::PBSpline{T}, i::PBSplineIndex, x::T, n::Int = b.n, nᵢ::Int = nnodes(b)) where {T}
    x *= nᵢ; x -= (i-1)
    # periodicity
    while x < 0  x += nᵢ end
    while x > nᵢ x -= nᵢ end
    return _unitpbspline(n,x)
end

# evaluation of "unit" spline with support (0,n+1)
function _unitpbspline(n::Int, x::T) where {T}
    if n == 0
        return zero(T) <= x < one(T) ? one(T) : zero(T)
    elseif x > n+1
        return zero(T)
    else
        return x/n * _unitpbspline(n-1,x) + (n+1-x)/n * _unitpbspline(n-1,x-1)
    end
end

function BasisFunctions.unsafe_eval_element(b::PBSpline, i, x)
    _pbspline(b, i, x)
end

function BasisFunctions.unsafe_eval_element_derivative(b::PBSpline, i, x, nᵢ::Int = nnodes(b))
    x *= nᵢ; x -= i-1
    # periodicity
    while x < 0  x += nᵢ end
    while x > nᵢ x -= nᵢ end
    nᵢ * (_unitpbspline(degree(b)-1, x) - _unitpbspline(degree(b)-1, x-1 < 0 ? x-1+nᵢ : x-1) )
end
