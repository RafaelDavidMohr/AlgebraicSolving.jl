# Benchmark for `equidimensional_decomposition` on the "Themos-net" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving

const NAME = "Themos-net"

function build_system(K)
    _, (x, y, a, b, c, d) = polynomial_ring(K, ["x", "y", "a", "b", "c", "d"])
    return [
        10000*d - 323536,
        x^2 - 96*x + y^2 - 178*y - 2304*a^2 - 8544*a*b - 2688*a*c - 7921*b^2 - 4984*b*c - 784*c^2 - d^2 + 11009,
        x^2 - 154*x + y^2 - 6*y - 5929*a^2 - 462*a*b - 5698*a*c - 9*b^2 - 222*b*c - 1369*c^2 - d^2 + 7307,
        x^2 - 98*x + y^2 - 46*y - 2401*a^2 - 2254*a*b - 5586*a*c - 529*b^2 - 2622*b*c - 3249*c^2 - d^2 + 6179,
        x*a + y*b,
        a^2 + b^2 + c^2 - 1
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:6])
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
