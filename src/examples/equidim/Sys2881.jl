# Benchmark for `equidimensional_decomposition` on the "Sys2881" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving

const NAME = "Sys2881"

function build_system(K)
    _, (x_1, x_10, x_11, x_12, x_13, x_14, x_15, x_2, x_3, x_4, x_5, x_6, x_7, x_8, x_9) = polynomial_ring(K, ["x_1", "x_10", "x_11", "x_12", "x_13", "x_14", "x_15", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9"])
    return [
        -x_13^2,
        -x_15^2,
        -x_3^2,
        -2*x_10*x_13,
        -2*x_12*x_15,
        -2*x_13*x_14,
        -2*x_14*x_15,
        -2*x_2*x_3,
        -2*x_3*x_6,
        -2*x_1*x_2 + x_5,
        -2*x_11*x_15 - 2*x_12*x_14,
        -2*x_10*x_14 - 2*x_11*x_13,
        -2*x_13*x_15 - x_14^2,
        -2*x_2*x_6 - 2*x_3*x_5,
        -x_1^2 + x_4 - 2,
        -2*x_1*x_4 + 2*x_1 + 2*x_7,
        -x_10^2 - x_13^2 - 2*x_13*x_7,
        -x_12^2 - x_15^2 - 2*x_15*x_9,
        -2*x_1*x_3 - x_2^2 + x_6,
        -x_3^2 - 2*x_3*x_9 - x_6^2,
        -2*x_10*x_15 - 2*x_11*x_14 - 2*x_12*x_13,
        -2*x_12*x_15 - 2*x_12*x_9 - 2*x_15*x_6,
        -2*x_12*x_3 - 2*x_3*x_6 - 2*x_6*x_9,
        -2*x_10*x_13 - 2*x_10*x_7 - 2*x_13*x_4 + 6*x_13,
        -2*x_2*x_3 - 2*x_2*x_9 - 2*x_3*x_8 - 2*x_5*x_6,
        -2*x_12*x_3 - 2*x_12*x_9 - 2*x_15*x_6 - 2*x_6*x_9,
        -2*x_1*x_5 - 2*x_2*x_4 + 5*x_2 + 2*x_8,
        -2*x_10*x_11 - 2*x_13*x_14 - 2*x_13*x_8 - 2*x_14*x_7,
        -2*x_11*x_12 - 2*x_14*x_15 - 2*x_14*x_9 - 2*x_15*x_8,
        -x_12^2 - 2*x_12*x_6 - 2*x_15*x_3 - 2*x_15*x_9 - x_9^2,
        -2*x_12*x_6 - 2*x_15*x_3 - 2*x_3*x_9 - x_6^2 - x_9^2,
        -2*x_1*x_6 - 2*x_2*x_5 - 2*x_3*x_4 + 8*x_3 + 2*x_9,
        -2*x_1*x_10 - 2*x_1*x_4 + 2*x_1 + 4*x_13 - 2*x_4*x_7 + 6*x_7,
        -2*x_1*x_13 - x_10^2 - 2*x_10*x_4 + 5*x_10 - 2*x_13*x_7 - x_7^2,
        -x_1^2 - 2*x_1*x_7 + 3*x_10 - x_4^2 + 4*x_4 - 4,
        -2*x_11*x_15 - 2*x_11*x_9 - 2*x_12*x_14 - 2*x_12*x_8 - 2*x_14*x_6 - 2*x_15*x_5,
        -2*x_1*x_10 - 2*x_10*x_7 - 2*x_13*x_4 + 10*x_13 - 2*x_4*x_7 + 4*x_7,
        -2*x_11*x_3 - 2*x_12*x_2 - 2*x_2*x_6 - 2*x_3*x_5 - 2*x_5*x_9 - 2*x_6*x_8,
        -2*x_1*x_2 - 2*x_1*x_8 + 3*x_11 - 2*x_2*x_7 - 2*x_4*x_5 + 7*x_5,
        -2*x_10*x_14 - 2*x_10*x_8 - 2*x_11*x_13 - 2*x_11*x_7 - 2*x_13*x_5 - 2*x_14*x_4 + 9*x_14,
        -2*x_10*x_12 - x_11^2 - 2*x_13*x_15 - 2*x_13*x_9 - x_14^2 - 2*x_14*x_8 - 2*x_15*x_7,
        -2*x_1*x_13 - 2*x_1*x_7 - 2*x_10*x_4 + 8*x_10 - x_4^2 + 3*x_4 - x_7^2 - 2,
        -2*x_11*x_3 - 2*x_11*x_9 - 2*x_12*x_2 - 2*x_12*x_8 - 2*x_14*x_6 - 2*x_15*x_5 - 2*x_5*x_9 - 2*x_6*x_8,
        -2*x_11*x_12 - 2*x_11*x_6 - 2*x_12*x_5 - 2*x_14*x_3 - 2*x_14*x_9 - 2*x_15*x_2 - 2*x_15*x_8 - 2*x_8*x_9,
        -2*x_11*x_6 - 2*x_12*x_5 - 2*x_14*x_3 - 2*x_15*x_2 - 2*x_2*x_9 - 2*x_3*x_8 - 2*x_5*x_6 - 2*x_8*x_9,
        -2*x_1*x_3 - 2*x_1*x_9 + 3*x_12 - x_2^2 - 2*x_2*x_8 - 2*x_3*x_7 - 2*x_4*x_6 - x_5^2 + 10*x_6,
        -2*x_1*x_11 - 2*x_1*x_5 - 2*x_10*x_2 + 4*x_14 - 2*x_2*x_4 + 5*x_2 - 2*x_4*x_8 - 2*x_5*x_7 + 9*x_8,
        -2*x_1*x_14 - 2*x_10*x_11 - 2*x_10*x_5 - 2*x_11*x_4 + 8*x_11 - 2*x_13*x_2 - 2*x_13*x_8 - 2*x_14*x_7 - 2*x_7*x_8,
        -2*x_10*x_15 - 2*x_10*x_9 - 2*x_11*x_14 - 2*x_11*x_8 - 2*x_12*x_13 - 2*x_12*x_7 - 2*x_13*x_6 - 2*x_14*x_5 - 2*x_15*x_4 + 12*x_15,
        -2*x_1*x_14 - 2*x_1*x_8 - 2*x_10*x_5 - 2*x_11*x_4 + 11*x_11 - 2*x_13*x_2 - 2*x_2*x_7 - 2*x_4*x_5 + 6*x_5 - 2*x_7*x_8,
        -2*x_1*x_11 - 2*x_10*x_2 - 2*x_10*x_8 - 2*x_11*x_7 - 2*x_13*x_5 - 2*x_14*x_4 + 13*x_14 - 2*x_4*x_8 - 2*x_5*x_7 + 7*x_8,
        -2*x_1*x_12 - 2*x_1*x_6 - 2*x_10*x_3 - 2*x_11*x_2 + 4*x_15 - 2*x_2*x_5 - 2*x_3*x_4 + 8*x_3 - 2*x_4*x_9 - 2*x_5*x_8 - 2*x_6*x_7 + 12*x_9,
        -2*x_1*x_15 - 2*x_10*x_12 - 2*x_10*x_6 - x_11^2 - 2*x_11*x_5 - 2*x_12*x_4 + 11*x_12 - 2*x_13*x_3 - 2*x_13*x_9 - 2*x_14*x_2 - 2*x_14*x_8 - 2*x_15*x_7 - 2*x_7*x_9 - x_8^2,
        -2*x_1*x_12 - 2*x_10*x_3 - 2*x_10*x_9 - 2*x_11*x_2 - 2*x_11*x_8 - 2*x_12*x_7 - 2*x_13*x_6 - 2*x_14*x_5 - 2*x_15*x_4 + 16*x_15 - 2*x_4*x_9 - 2*x_5*x_8 - 2*x_6*x_7 + 10*x_9,
        -2*x_1*x_15 - 2*x_1*x_9 - 2*x_10*x_6 - 2*x_11*x_5 - 2*x_12*x_4 + 14*x_12 - 2*x_13*x_3 - 2*x_14*x_2 - 2*x_2*x_8 - 2*x_3*x_7 - 2*x_4*x_6 - x_5^2 + 9*x_6 - 2*x_7*x_9 - x_8^2
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:15])
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
