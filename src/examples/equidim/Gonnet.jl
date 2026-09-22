# Benchmark for `equidimensional_decomposition` on the "Gonnet" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving

const NAME = "Gonnet"

function build_system(K)
    _, (x, y, z, t, u, v, w, a, b, c, d, e, f, g, h, i, j) = polynomial_ring(K, ["x", "y", "z", "t", "u", "v", "w", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j"])
    return [
        a*g,
        w + a*h + b*g + d*f,
        c*i,
        w*f,
        t + a*i + c*g + c*j + c + d*i + e*i,
        y + u + v + w*i + a*g + a*j + a + b*i + c*f + c*h + d*g + d*j + d + e*g + e*j + e,
        x + z + w*g + w*j + w + a*f + a*h + b*g + b*j + b + d*f + d*h + e*f + e*h - 1,
        w*f + w*h + b*f + b*h,
        2*w*f + w*h + b*f,
        x + 2*w*g + w*j + w + 2*a*f + a*h + b*g + d*f + e*f,
        y + w*i + 2*a*g + a*j + a + c*f + d*g + e*g,
        a*i + c*g,
        w*g + a*f,
        w*h + b*f + 2*b*h,
        z + w + a*h + b*g + b*j + 2*b + d*f + 2*d*h + e*h,
        u + a + b*i + c*h + d*g + d*j + 2*d + e,
        c + d*i,
        w*h + b*f,
        a + d*g
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:17])
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
