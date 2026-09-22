# Benchmark for `equidimensional_decomposition` on the "Sys2943" system.

using Pkg
Pkg.activate(joinpath(@__DIR__, "..", "..", ".."))

using AlgebraicSolving
using FileWatching.Pidfile: mkpidlock

const NAME = "Sys2943"

function build_system(K)
    _, (Y, c12, c13, c23, d12, d13, d23, s12, s13, s23, t12, t13, t23) = polynomial_ring(K, ["Y", "c12", "c13", "c23", "d12", "d13", "d23", "s12", "s13", "s23", "t12", "t13", "t23"])
    return [
        c12*t12 - s12,
        c13*t13 - s13,
        c23*t23 - s23,
        c12^2 + s12^2 - 1,
        c13^2 + s13^2 - 1,
        c23^2 + s23^2 - 1,
        -2*c12*d13*d23 - d12^2 + d13^2 + d23^2,
        -2*c13*d12*d23 + d12^2 - d13^2 + d23^2,
        -2*c23*d12*d13 + d12^2 + d13^2 - d23^2,
        Y*d12*d13*d23 - 1
    ]
end

# Decompose a tiny system first so that the timings below measure the
# computation rather than compilation. The solver specialises on the number of
# variables, so the warm up uses as many of them as the benchmark itself.
function warmup(K)
    _, w = polynomial_ring(K, ["w$(i)" for i in 1:13])
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
time_qq = @elapsed equidimensional_decomposition(Ideal(build_system(QQ)), info_level = 1)
if time_qq < 120
    println("### $(NAME): second run over QQ ###")
    time_qq = @elapsed equidimensional_decomposition(Ideal(build_system(QQ)), info_level = 1)
end

println("### $(NAME): equidimensional_decomposition over GF($(prime)) ###")
time_gf = @elapsed equidimensional_decomposition(Ideal(build_system(GF(prime))), info_level = 1)
if time_gf < 120
    println("### $(NAME): second run over GF($(prime)) ###")
    time_gf = @elapsed equidimensional_decomposition(Ideal(build_system(GF(prime))), info_level = 1)
end

ratio = time_qq / time_gf

println()
println("### $(NAME) timings ###")
println("QQ:          $(time_qq) s")
println("GF($(prime)): $(time_gf) s")
println("QQ / GF:     $(ratio)")

# Collect the timings of every benchmark in this folder in one table. Rows are
# appended, so re-running a benchmark adds a row rather than replacing one.
# Concurrent runs would otherwise lose and interleave rows, hence the lock.
csv = joinpath(@__DIR__, "timings.csv")
mkpidlock(csv * ".lock") do
    write_header = !isfile(csv)
    open(csv, "a") do io
        write_header && println(io, "system,nvars,neqns,time_qq,time_gf,ratio_qq_gf")
        println(io, join((NAME, 13, 10,
                          round(time_qq, digits = 6),
                          round(time_gf, digits = 6),
                          round(ratio, digits = 4)), ","))
    end
end
println("appended to $(csv)")
