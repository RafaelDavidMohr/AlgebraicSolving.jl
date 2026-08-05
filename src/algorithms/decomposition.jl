@doc Markdown.doc"""
    equations(X::LocallyClosedSet)

Given a locally closed set `X` of the form $V(F) \ V(g_1 \cdot \dots \cdot g_r)$,
return the list of polynomials $F$.
"""
equations(X::LocallyClosedSet) = X.eqns

@doc Markdown.doc"""
    inequations(X::LocallyClosedSet)

Given a locally closed set `X` of the form $V(F) \ V(g_1 \cdot
\dots \cdot g_r)$, return the list of polynomials $g_1, \dots g_r$.
"""
inequations(X::LocallyClosedSet) = X.ineqns

@doc Markdown.doc"""
    groebner_basis(X::LocallyClosedSet)

Return a Gröbner basis for an ideal whose zeros coincide with the
Zariski closure of `X`. If no Gröbner basis is present, one is computed
from scratch.
"""
function groebner_basis(X::LocallyClosedSet)
    !isempty(X.gb) && return X.gb
    gb = saturate(X.eqns, X.ineqns)
    X.gb = gb
    return gb
end

@doc Markdown.doc"""
    equidimensional_decomposition(I::Ideal{T}, info_level::Int=0) where {T <: MPolyRingElem}

Given a polynomial ideal `I`, return a list of locally closed sets
`dec` s.t.  each locally closed set in `dec` is equidimensional
(i.e. its Zariski closure has irreducible components of only one
dimension) and s.t. the algebraic set defined by `I` equals the union of
the locally closed sets in `dec`.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `I::Ideal{T} where T <: MpolyElem`: input ideal.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)

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

    F = I.gens
    Fhom = homogenize(F)
    sort!(Fhom, by = p -> total_degree(p))
    r = ModularRegistry(T[])
    cells = LocClosedSet{FqMPolyRingElem}[]
    while !is_finished(r)
        cells = _sig_decomp(Fhom, r, info_level = info_level)
    end
    res = LocallyClosedSet{T}[]
    R = parent(I)
    for cell in cells
        append!(res, get_output_cells(cell, R, F, r))
    end
    return res
end
