
using FastGaussQuadrature

get_nodes(::Val{:gauss_legendre}, n, T=Float64) = convert(Array{T}, gausslegendre(n)[1])
get_nodes(::Val{:lobatto_legendre}, n, T=Float64) = convert(Array{T}, gausslobatto(n)[1])

function get_nodes(::Val{:gauss_chebyshev}, n, T=Float64; k=1)
    x = zeros(T, n)
    if k == 1
        for i in eachindex(x)
            x[i] = sin( π*(n-2i+1) / (2n) )
        end
        return x[end:-1:1]
    elseif k == 2
        for i in eachindex(x)
            x[i] = cos( π*i / (n-1) )
        end
        return x
    end
end

function get_nodes(::Val{:lobatto_chebyshev}, n, T=Float64)
    x = zeros(T, n)
    for i in eachindex(x)
        x[i] = cos( (i-1)*π / (n-1) )
    end
    return x[end:-1:1]
end
