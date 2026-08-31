@doc Markdown.doc"""
    equidimensional_decomposition(I::Ideal{T}, info_level::Int=0) where {T <: MPolyRingElem}

Given a polynomial ideal `I`, return a list of locally closed sets
`dec` s.t.  each locally closed set in `dec` is equidimensional
(i.e. its Zariski closure has irreducible components of only one
dimension) and s.t. the algebraic set defined by `I` equals the union of
the locally closed sets in `dec`.

When the ground field of `I` is the rational numbers, a multi-modular method is used.

# Arguments
- `I::Ideal{T} where T <: MPolyRingElem`: input ideal.
- `info_level::Int=0`: info level printout: off (`0`, default), details about multi-modular computation (`1`), computational details (`2`)

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(GF(65521), ["x", "y", "z"])
(Multivariate polynomial ring in 3 variables over GF(65521), FqMPolyRingElem[x, y, z])

julia> I = Ideal([x*y - x*z, x*z^2 - x*z, x^2*z - x*z])
FqMPolyRingElem[x*y + 65520*x*z, x*z^2 + 65520*x*z, x^2*z + 65520*x*z]

julia> equidimensional_decomposition(I)
3-element Vector{AlgebraicSolving.LocallyClosedSet{FqMPolyRingElem}}:
 V(x*y + 65520*x*z) \ V(y + 65520*z)
 V(y + 65520*z, x*z^2 + 65520*x*z) \ V(x * (z + 65520))
 V(y + 65520*z, x*z^2 + 65520*x*z, x^2*z + 65520*x*z, z + 65520) \ V(x * z)
```
"""
function equidimensional_decomposition(I::Ideal{T};
                                       info_level::Int=0) where {T <: MPolyRingElem}
    return _equidimensional_decomposition(I, info_level = info_level)
end

function _equidimensional_decomposition(I::Ideal{T};
                                        info_level::Int=0) where {T <: FqMPolyRingElem}

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    r = ModularRegistry(T[])
    cells = _sig_decomp(Fhom, r, info_level = info_level - 1)
    res = LocallyClosedSet{T}[]
    R = parent(I)
    for cell in cells
        append!(res, get_output_cells(cell, R, r))
    end
    return res
end

function _equidimensional_decomposition(I::Ideal{T};
                                        info_level::Int=0) where {T <: QQMPolyRingElem}

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    Rhom = parent(first(Fhom))
    r = ReconstructRegistry(Rhom, ReconstructPol[], 1,
                            Int32[], Int32(0))
    cells = LocClosedSet{FqMPolyRingElem}[]
    cnt = 0
    while !is_finished(r)
        cnt += 1
        p = Int32(rand_bits_prime(ZZ, 31))
            
        new_prime!(r, p)
        S, _ = polynomial_ring(GF(p), ["x$i" for i in 1:ngens(Rhom)],
                               internal_ordering = :degrevlex)
        Fhomp = [reduce_mod_p(f, S) for f in Fhom]
        cells = _sig_decomp(Fhomp, r, info_level = info_level - 1)
        if info_level >= 1
            npols = length(r.pols)
            nstable = length(findall(p -> p.is_stable, r.pols))
            @info "decomposition $cnt with prime $p, $nstable / $npols finished"
        end
    end
    res = LocallyClosedSet{T}[]
    R = parent(I)
    for cell in cells
        append!(res, get_output_cells(cell, R, Fhom, r))
    end
    info_level >= 1 && @info "needed $cnt primes"
    return res
end

