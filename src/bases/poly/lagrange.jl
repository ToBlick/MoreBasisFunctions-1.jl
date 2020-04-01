#######################
# The Lagrange basis
#######################

using Polynomials: Poly, polyint
using BasisFunctions: ChebyshevInterval, PolynomialBasis
using BasisFunctions: hasderivative, hasantiderivative, ordering, support

const LagrangeInterval = ChebyshevInterval
const LagrangeIndex = NativeIndex{:lagrange}

"""
A basis of the Lagrange polynomials `l_i(x) = ∏_(j,i≠j) (x - ξ^j) / (ξ^i - ξ^j)`
on the interval [-1,+1].
"""
struct Lagrange{S,T} <: PolynomialBasis{T}
    n :: Int
    nodes :: ScatteredGrid{T}

    denom::Vector{T}    # denominator of Lagrange polynomials
    diffs::Matrix{T}    # inverse of differences between nodes
    vdminv::Matrix{T}   # inverse Vandermonde matrix

    function Lagrange{S,T}(nodes::ScatteredGrid{T}) where {S,T}
        local p::T

        local ξ = nodes.points
        local n = length(ξ)

        @assert minimum(ξ) ≥ leftendpoint(support(nodes))
        @assert maximum(ξ) ≤ rightendpoint(support(nodes))

        denom = ones(n)
        diffs = zeros(n,n)

        for i in eachindex(ξ)
            p = 1
            for j in eachindex(ξ)
                diffs[i,j] = 1 / (ξ[i] - ξ[j])
                if i ≠ j
                    denom[i] *= diffs[i,j]
                end
            end
        end

        new(n, nodes, denom, diffs, vandermonde_matrix_inverse(ξ))
    end
end

Lagrange{S}(nodes::ScatteredGrid{T}) where {S,T} = Lagrange{S,T}(nodes)
Lagrange(nodes::ScatteredGrid) = Lagrange{:custom}(nodes)
Lagrange{S}(ξ::Vector) where {S} = Lagrange{S}(ScatteredGrid(ξ, LagrangeInterval()))
Lagrange(ξ::Vector) = Lagrange(ScatteredGrid(ξ, LagrangeInterval()))
Lagrange{S}(n::Int) where {S} = Lagrange{S}(get_nodes(Val(S), n))
Lagrange{S,T}(n::Int) where {S,T} = Lagrange{S}(get_nodes(Val(S), n, T))

const LagrangeGaussLegendre = Lagrange{:gauss_legendre}
const LagrangeLobattoLegendre = Lagrange{:lobatto_legendre}
const LagrangeGaussChebyshev = Lagrange{:gauss_chebyshev}
const LagrangeLobattoChebyshev = Lagrange{:lobatto_chebyshev}

# Convenience constructor: map the Lagrange basis to the interval [a,b]
Lagrange{S}(n::Int, a::Number, b::Number) where {S} = rescale(Lagrange{S}(n), a, b)


nodes(b::Lagrange)  = b.nodes.points
nnodes(b::Lagrange) = b.n
degree(b::Lagrange) = b.n-1

BasisFunctions.native_index(b::Lagrange, idx) = LagrangeIndex(idx)
BasisFunctions.linear_index(b::Lagrange, idx) = LagrangeIndex(idx)
BasisFunctions.ordering(b::Lagrange) = Base.OneTo(nnodes(b))
BasisFunctions.support(b::Lagrange) = support(b.nodes)

BasisFunctions.interpolation_grid(b::Lagrange{S,T}) where {S,T} = ScatteredGrid(get_nodes(Val(S), length(b), T), LagrangeInterval())

BasisFunctions.similar(::Lagrange{S}, ::Type{T}, n::Int) where {S,T} = Lagrange{S,T}(n)

Base.size(b::Lagrange) = (b.n,)


function BasisFunctions.unsafe_eval_element(b::Lagrange{S,T}, idx::LagrangeIndex, x::T) where {S,T}
    local y::T = 1
    local ξ = nodes(b)
    for i in 1:length(ξ)
        i ≠ idx ? y *= x - ξ[i] : nothing
    end
    y * b.denom[value(idx)]
end


function BasisFunctions.unsafe_eval_element_derivative(b::Lagrange{S,T}, idx::LagrangeIndex, x::T) where {S,T}
    local y::T = 0
    local z::T
    local ξ = nodes(b)
    local d = b.diffs

    for l in 1:nnodes(b)
        if l ≠ idx
            z = d[idx,l]
            for i in 1:nnodes(b)
                i ≠ idx && i ≠ l ? z *= (x - ξ[i]) * d[idx,i] : nothing
            end
            y += z
        end
    end
    y
end


function unsafe_eval_element_antiderivative(b::Lagrange{S,T}, idx::LagrangeIndex, x::T) where {S,T}
    local y = zero(nodes(b))
    y[value(idx)] = 1
    lint = polyint(Poly(b.vdminv*y))
    return lint(x) - lint(leftendpoint(support(b)))
end
