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
    p :: Int                    # spline degree
    nodes :: ScatteredGrid{T}   # knot points
    # n :: Int                    # number of splines
    # L :: T                      # domain length
    # h :: T                      # L/n

    function PBSpline(p::Int, nodes::ScatteredGrid{T}) where {T}
        @assert length(nodes) > p
        @assert check_equispaced(nodes.points)

        new{T}(p,nodes)
    end
end

# constructors for unit interval grid
PBSpline(p::Int, n::Int) = PBSpline(p, get_pbspline_nodes(Float64,n))
PBSpline(p::Int, ξ::Vector{T}) where {T} = PBSpline(p, ScatteredGrid(ξ, PBSplineInterval{T}()))

# n splines of degree p=n-1
PBSpline(n::Int) = PBSpline(n-1, get_pbspline_nodes(Float64,n))
PBSpline{T}(n::Int) where {T} = PBSpline(n-1, get_pbspline_nodes(T,n))
PBSpline(ξ::Vector{T}) where {T} = PBSpline(length(ξ)-1, ξ)

# constructor for interval [a,b)
PBSpline(p, ξ, a::Number, b::Number) = rescale(PBSpline(p, ξ), a, b)

# basis properties
nodes(b::PBSpline)  = b.nodes.points
nnodes(b::PBSpline) = length(b.nodes.points)
degree(b::PBSpline) = b.p
Base.size(b::PBSpline) = (nnodes(b),)

# convenience
L_domain(b::PBSpline) = nodes(b)[end] - nodes(b)[1]
invh_elements(b::PBSpline) = nnodes(b)/L_domain(b)

BasisFunctions.native_index(b::PBSpline, idx) = PBSplineIndex(idx)
BasisFunctions.linear_index(b::PBSpline, idx) = PBSplineIndex(idx)
BasisFunctions.support(b::PBSpline{T}) where {T} = PBSplineInterval{T}()

BasisFunctions.interpolation_grid(b::PBSpline{T}) where {T} = ScatteredGrid(get_pbspline_nodes(T,length(b)), PBSplineInterval{T}())

BasisFunctions.similar(::PBSpline, ::Type{T}, n::Int) where {T} = PBSpline{T}(n)
BasisFunctions.similar(::PBSpline, ::Type{T}, ξ::Vector{T}) where {T} = PBSpline(ξ)

# evaluation at point x
function _pbspline(b::PBSpline{T}, i::PBSplineIndex, x::T,
                    n::Int = nnodes(b), invh::T = invh_elements(b)) where {T}
    x -= nodes(b)[1]; x *= invh; x -= (i-1)
    # periodicity
    while x < zero(T)  x += n end
    while x > n        x -= n end
    return _unitpbspline(b.p,x)
end

# evaluation of "unit" spline with support (0,n+1)
function _unitpbspline(p::Int, x::T) where {T}
    if p == 0
        return zero(T) <= x < one(T) ? one(T) : zero(T)
    elseif x > p+1
        return zero(T)
    else
        return x/p * _unitpbspline(p-1,x) + (p+1-x)/p * _unitpbspline(p-1,x-1)
    end
end

function BasisFunctions.unsafe_eval_element(b::PBSpline, i, x)
    _pbspline(b, i, x)
end

function _pbspline_derivative(b::PBSpline{T}, i::PBSplineIndex, x::T,
                    n::Int = nnodes(b), invh::T = invh_elements(b)) where {T}
    x -= nodes(b)[1]; x *= invh; x -= (i-1)
    # periodicity
    while x < zero(T)  x += n end
    while x > n        x -= n end
    invh * (_unitpbspline(b.p-1, x) - _unitpbspline(b.p-1, x-1 < 0 ? x-1+n : x-1) )
end

function BasisFunctions.unsafe_eval_element_derivative(b::PBSpline, i, x, nᵢ::Int = nnodes(b))
    _pbspline_derivative(b, i, x)
end
