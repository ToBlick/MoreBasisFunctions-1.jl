#######################
# The Bernstein basis
#######################

using OffsetArrays

using BasisFunctions: UnitInterval, PolynomialBasis
using BasisFunctions: hasderivative, hasantiderivative, support

const BernsteinInterval = UnitInterval

struct Bernstein{T} <: PolynomialBasis{T,T}
    n :: Int
    ξ :: OffsetVector{T,Vector{T}}

    function Bernstein{T}(ξ) where {T}
        new(length(ξ), OffsetVector(ξ, 0:length(ξ)-1))
    end
end

function Bernstein(ξ::Vector{T}) where {T}
    Bernstein{T}(ξ)
end

# Convenience constructor: map the Bernstein basis to the interval [a,b]
Bernstein(ξ, a::Number, b::Number) = rescale(Bernstein(ξ), a, b)


nodes(b::Bernstein)  = b.ξ.parent
nnodes(b::Bernstein) = b.n
degree(b::Bernstein) = b.n-1

BasisFunctions.native_index(b::Bernstein, idxn) = PolynomIndex(idxn)
BasisFunctions.hasderivative(b::Bernstein) = true
BasisFunctions.hasantiderivative(b::Bernstein) = false
BasisFunctions.support(b::Bernstein) = BernsteinInterval()

Base.size(b::Bernstein) = b.n


function _bernstein(b::Bernstein{T}, i::Int, n::Int, x::T) where {T}
    if i < 0 || i > n
        return zero(T)
    else
        if n == 0
            return one(T)
        else
            return _bernstein(b, i, n-1, x) * (1-x) + _bernstein(b, i-1, n-1, x) * x
        end
    end
end


function BasisFunctions.unsafe_eval_element(b::Bernstein{T}, i::PolynomIndex, x::T) where {T}
    @assert i ≥ 1 && i ≤ nnodes(b)
    _bernstein(b, i-1, nnodes(b)-1, x)
end


function BasisFunctions.unsafe_eval_element_derivative(b::Bernstein{T}, i::PolynomIndex, x::T) where {T}
    @assert i ≥ 1 && i ≤ nnodes(b)
    (nnodes(b)-1) * ( _bernstein(b, i-2, nnodes(b)-2, x) - _bernstein(b, i-1, nnodes(b)-2, x) )
end
