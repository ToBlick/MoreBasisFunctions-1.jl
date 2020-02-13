
using BasisFunctions
using MoreBasisFunctions
using DomainSets
using FrameFun


f(x) = cos(π*x)

L = LagrangeGaussLegendre(21)
FL = Fun(f, L, support(L))
println(FL(0.0))

B = Bernstein(21)
FB = Fun(f, B, support(B))
println(FB(0.5))
