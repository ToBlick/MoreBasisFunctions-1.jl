
using MoreBasisFunctions: vandermonde_matrix, vandermonde_matrix_inverse

n = 11
x = collect(range(0., stop=1., length=n))
y = rand(n)

V = vandermonde_matrix(x)
W = vandermonde_matrix_inverse(x)

x = V\y
a = W*y

@test maximum(abs.((a .- x) ./ a)) < eps(Float32)
