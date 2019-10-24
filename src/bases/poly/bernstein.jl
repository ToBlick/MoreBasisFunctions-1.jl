#######################
# The Bernstein basis
#######################

using OffsetArrays

using BasisFunctions: UnitInterval, PolynomialBasis
using BasisFunctions: hasderivative, hasantiderivative, support

const BernsteinInterval = UnitInterval
const BernsteinIndex = NativeIndex{:bernstein}

get_bernstein_nodes(T,n) = collect(LinRange{T}(0.0, 1.0, n))

struct Bernstein{T} <: PolynomialBasis{T,T}
    n :: Int
    nodes :: ScatteredGrid{T}

    function Bernstein(nodes::ScatteredGrid{T}) where {T}
        new{T}(length(nodes.points), nodes)
    end
    function Bernstein(ξ::Vector{T}) where {T}
        Bernstein(ScatteredGrid(ξ, BernsteinInterval{T}()))
    end
end

Bernstein(n::Int) = Bernstein(get_bernstein_nodes(Float64,n))
Bernstein{T}(n::Int) where {T} = Bernstein(get_bernstein_nodes(T,n))

# Convenience constructor: map the Bernstein basis to the interval [a,b]
Bernstein(ξ, a::Number, b::Number) = rescale(Bernstein(ξ), a, b)


nodes(b::Bernstein)  = b.nodes.points
nnodes(b::Bernstein) = b.n
degree(b::Bernstein) = b.n-1

BasisFunctions.native_index(b::Bernstein, idx) = BernsteinIndex(idx)
BasisFunctions.linear_index(b::Bernstein, idx) = BernsteinIndex(idx)
BasisFunctions.support(b::Bernstein{T}) where {T} = BernsteinInterval{T}()

BasisFunctions.interpolation_grid(b::Bernstein{T}) where {T} = ScatteredGrid(get_bernstein_nodes(T,length(b)), BernsteinInterval{T}())

BasisFunctions.similar(::Bernstein, ::Type{T}, n::Int) where {T} = Bernstein{T}(n)
BasisFunctions.similar(::Bernstein, ::Type{T}, ξ::Vector{T}) where {T} = Bernstein(ξ)

Base.size(b::Bernstein) = (b.n,)


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
    _bernstein(b, i-1, nnodes(b)-1, x)
end


function BasisFunctions.unsafe_eval_element_derivative(b::Bernstein, i, x)
    (nnodes(b)-1) * ( _bernstein(b, i-2, nnodes(b)-2, x) - _bernstein(b, i-1, nnodes(b)-2, x) )
end
