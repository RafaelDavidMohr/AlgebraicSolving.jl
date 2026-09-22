# Benchmark for `equidimensional_decomposition` on the "8-3-config-Li" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving

const NAME = "8-3-config-Li"

function build_system(K)
    _, (x12, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1) = polynomial_ring(K, ["x12", "x11", "x10", "x9", "x8", "x7", "x6", "x5", "x4", "x3", "x2", "x1"])
    return [
        x12*x3 - x8*x7,
        x11*x5 - x10*x6,
        x9*x3 - x9*x1 - x8*x4 + x8*x1,
        x12*x6 - x12*x1 - x11*x7 + x11*x1,
        x10*x3 - x10*x2 - x8*x5 + x8*x2,
        x11*x4 - x11*x2 - x9*x6 + x9*x2,
        x12*x5 - x12*x4 - x10*x7 + x10*x4 + x9*x7 - x9*x5
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:12])
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
