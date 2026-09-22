# Benchmark for `equidimensional_decomposition` on the "Leykin-1" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving

const NAME = "Leykin-1"

function build_system(K)
    _, (x4, x3, x2, x1, a4, a3, a2, a1) = polynomial_ring(K, ["x4", "x3", "x2", "x1", "a4", "a3", "a2", "a1"])
    return [
        -7*x4*a4 - 3*x3*a3 + x2*a2 + 5*x1*a1,
        -4*x4*x2*a4 - 5*x3^2*a4 - 6*x3*x2*a3 - 5*x3*x1*a2 + 2*x2^2*a2,
        -2*x4*x1*a4 + 5*x3*x2*a4 + 2*x3*x1*a3 + x2*x1*a2,
        x3*x1*a4 + 2*x2^2*a4 + 2*x2*x1*a3 + x1^2*a2,
        7*x4*x3*a4 + 10*x4*x2*a3 + 5*x4*x1*a2 - 2*x3^2*a3 + 4*x3*x2*a2 + 10*x2*x1*a1,
        x4*x3*x1^2 - x4*x2^2*x1 - x3^2*x2*x1 + x3*x2^3
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:8])
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
