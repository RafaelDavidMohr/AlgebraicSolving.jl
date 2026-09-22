# Benchmark for `equidimensional_decomposition` on the "mohab2" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving

const NAME = "mohab2"

function build_system(K)
    _, (u, v, w) = polynomial_ring(K, ["u", "v", "w"])
    return [
        160*u^4*v^3*w^3 - 64*u^3*v^4*w^4 - 288*u^3*v^4*w^2 + 108*u^3*v^4 - 288*u^3*v^2*w^4 - 24*u^3*v^2*w^2 + 108*u^3*w^4 + 96*u^2*v^5*w^3 + 96*u^2*v^3*w^5 + 552*u^2*v^3*w^3 - 216*u^2*v^3*w - 216*u^2*v*w^3 + 72*u^2*v*w - 144*u*v^4*w^4 - 12*u*v^4*w^2 - 12*u*v^2*w^4 + 168*u*v^2*w^2 - 36*u*v^2 - 36*u*w^2 + 8*u - 72*v^3*w^3 + 24*v^3*w + 24*v*w^3 - 8*v*w,
        96*u^5*v^2*w^3 - 64*u^4*v^3*w^4 - 288*u^4*v^3*w^2 + 108*u^4*v^3 - 144*u^4*v*w^4 - 12*u^4*v*w^2 + 160*u^3*v^4*w^3 + 96*u^3*v^2*w^5 + 552*u^3*v^2*w^3 - 216*u^3*v^2*w - 72*u^3*w^3 + 24*u^3*w - 288*u^2*v^3*w^4 - 24*u^2*v^3*w^2 - 12*u^2*v*w^4 + 168*u^2*v*w^2 - 36*u^2*v - 216*u*v^2*w^3 + 72*u*v^2*w + 24*u*w^3 - 8*u*w + 108*v^3*w^4 - 36*v*w^2 + 8*v,
        96*u^5*v^3*w^2 - 64*u^4*v^4*w^3 - 144*u^4*v^4*w - 288*u^4*v^2*w^3 - 12*u^4*v^2*w + 108*u^4*w^3 + 96*u^3*v^5*w^2 + 160*u^3*v^3*w^4 + 552*u^3*v^3*w^2 - 72*u^3*v^3 - 216*u^3*v*w^2 + 24*u^3*v - 288*u^2*v^4*w^3 - 12*u^2*v^4*w - 24*u^2*v^2*w^3 + 168*u^2*v^2*w - 36*u^2*w - 216*u*v^3*w^2 + 24*u*v^3 + 72*u*v*w^2 - 8*u*v + 108*v^4*w^3 - 36*v^2*w + 8*w
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:3])
    equidimensional_decomposition(Ideal([w[1]*w[2], w[1]*w[3], w[2]*w[3]]))
    return nothing
end

prime = Int32(AlgebraicSolving.Nemo.rand_bits_prime(ZZ, 31))
warmup(QQ)
warmup(GF(prime))

# A short first run can still be dominated by one off compilation that the warm
# up does not cover, so measure a second one and keep that timing. Past this
# threshold the computation dwarfs the overhead and the first timing stands.
println("### $(NAME): equidimensional_decomposition over QQ ###")
time_qq = @elapsed equidimensional_decomposition(Ideal(build_system(QQ)), info_level = 2)
if time_qq < 120
    println("### $(NAME): second run over QQ ###")
    time_qq = @elapsed equidimensional_decomposition(Ideal(build_system(QQ)), info_level = 2)
end

println("### $(NAME): equidimensional_decomposition over GF($(prime)) ###")
time_gf = @elapsed equidimensional_decomposition(Ideal(build_system(GF(prime))), info_level = 2)
if time_gf < 120
    println("### $(NAME): second run over GF($(prime)) ###")
    time_gf = @elapsed equidimensional_decomposition(Ideal(build_system(GF(prime))), info_level = 2)
end

println()
println("### $(NAME) timings ###")
println("QQ:          $(time_qq) s")
println("GF($(prime)): $(time_gf) s")
