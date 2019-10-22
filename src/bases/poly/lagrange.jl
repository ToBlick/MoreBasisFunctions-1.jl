#######################
# The Lagrange basis
#######################

using Polynomials: Poly, polyint
using BasisFunctions: ChebyshevInterval
using BasisFunctions: hasderivative, hasantiderivative, support

const LagrangeInterval = ChebyshevInterval

"""
A basis of the Lagrange polynomials `l_i(x) = ∏_(j,i≠j) (x - ξ^j) / (ξ^i - ξ^j)`
on the interval [-1,+1].
"""
struct Lagrange{T} <: PolynomialBasis{T,T}
    n :: Int
    grid :: ScatteredGrid{T}

    denom::Vector{T}
    diffs::Matrix{T}
    vdminv::Matrix{T}

    function Lagrange{T}(g) where {T}
        local p::T

        local ξ = g.points
        local n = length(ξ)

        @assert minimum(ξ) ≥ leftendpoint(support(g))
        @assert maximum(ξ) ≤ rightendpoint(support(g))

        denom = zeros(n)
        diffs = zeros(n,n)

        for i in 1:length(ξ)
            p = 1
            for j in 1:length(ξ)
                diffs[i,j] = ξ[i] - ξ[j]
                if i ≠ j
                    p *= diffs[i,j]
                end
            end
            denom[i] = 1/p
        end

        new(n, g, denom, diffs, vandermonde_matrix_inverse(ξ))
    end
end

function Lagrange(g::ScatteredGrid{T}) where {T}
    Lagrange{T}(g)
end

function Lagrange(ξ::Vector{T}) where {T}
    Lagrange(ScatteredGrid(ξ, LagrangeInterval()))
end

nodes(b::Lagrange)  = b.grid.points
nnodes(b::Lagrange) = b.n
degree(b::Lagrange) = b.n-1

BasisFunctions.native_index(b::Lagrange, idxn) = PolynomIndex(idxn)
BasisFunctions.hasderivative(b::Lagrange) = true
BasisFunctions.hasantiderivative(b::Lagrange) = true
BasisFunctions.support(b::Lagrange) = support(b.grid)

Base.size(b::Lagrange) = b.n


function BasisFunctions.unsafe_eval_element(b::Lagrange{T}, idx::PolynomIndex, x) where {T}
    local y::T = 1
    local ξ = b.grid.points
    for i in 1:length(ξ)
        i ≠ idx ? y *= x - ξ[i] : nothing
    end
    y * b.denom[idx]
end


function BasisFunctions.unsafe_eval_element_derivative(b::Lagrange{T}, idx::PolynomIndex, x) where {T}
    local y::T = 0
    local z::T
    local ξ = b.grid.points
    local d = b.diffs

    for l in 1:nnodes(b)
        if l ≠ idx
            z = 1 / d[idx,l]
            for i in 1:nnodes(b)
                i ≠ idx && i ≠ l ? z *= (x - ξ[i]) / d[idx,i] : nothing
            end
            y += z
        end
    end
    y
end


function unsafe_eval_element_antiderivative(b::Lagrange{T}, idx::PolynomIndex, x) where {T}
    local y = zero(b.grid.points)
    y[idx] = 1
    lint = polyint(Poly(b.vdminv*y))
    return lint(x) - lint(leftendpoint(support(b)))
end


similar(b::Lagrange, ::Type{T}, g::ScatteredGrid{T}) where {T} = Lagrange{T}(g)
