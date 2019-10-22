#######################
# The Bernstein basis
#######################

using OffsetArrays

using BasisFunctions: UnitInterval, PolynomialBasis
using BasisFunctions: hasderivative, hasantiderivative, support

const BernsteinInterval = UnitInterval
const BernsteinIndex = NativeIndex{:bernstein}

struct Bernstein{T} <: PolynomialBasis{T,T}
    n :: Int
    ξ :: OffsetVector{T,Vector{T}}

    function Bernstein(ξ::Vector{T}) where {T}
        new{T}(length(ξ), OffsetVector(ξ, 0:length(ξ)-1))
    end
end

Bernstein(n::Int) = Bernstein(collect(LinRange(0.0, 1.0, n)))
Bernstein{T}(n::Int) where {T} = Bernstein(collect(LinRange{T}(0.0, 1.0, n)))

# Convenience constructor: map the Bernstein basis to the interval [a,b]
Bernstein(ξ, a::Number, b::Number) = rescale(Bernstein(ξ), a, b)


nodes(b::Bernstein)  = b.ξ.parent
nnodes(b::Bernstein) = b.n
degree(b::Bernstein) = b.n-1

BasisFunctions.native_index(b::Bernstein, idxn) = BernsteinIndex(idxn)
BasisFunctions.hasderivative(b::Bernstein) = true
BasisFunctions.hasantiderivative(b::Bernstein) = false
BasisFunctions.support(b::Bernstein) = BernsteinInterval()

BasisFunctions.similar(::Bernstein, ::Type{T}, n::Int) where {T} = Bernstein{T}(n)
BasisFunctions.similar(::Bernstein, ::Type{T}, ξ::Vector{T}) where {T} = Bernstein(ξ)

Base.size(b::Bernstein) = b.n


function _bernstein(b::Bernstein{T}, i::BernsteinIndex, n::BernsteinIndex, x::T) where {T}
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

_bernstein(b::Bernstein, i, n, x) = _bernstein(b, native_index(b,i), native_index(b,n), x)


function BasisFunctions.unsafe_eval_element(b::Bernstein, i, x)
    @assert i ≥ 1 && i ≤ nnodes(b)
    _bernstein(b, i-1, nnodes(b)-1, x)
end


function BasisFunctions.unsafe_eval_element_derivative(b::Bernstein, i, x)
    @assert i ≥ 1 && i ≤ nnodes(b)
    (nnodes(b)-1) * ( _bernstein(b, i-2, nnodes(b)-2, x) - _bernstein(b, i-1, nnodes(b)-2, x) )
end
