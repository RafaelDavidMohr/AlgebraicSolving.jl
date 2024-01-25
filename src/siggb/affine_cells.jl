# for convenience
function affine_cell(F::Vector{T},
                     H::Vector{T}) where {T <: MPolyRingElem}
    return LocClosedSet{T}(F,H)
end

# for displaying locally closed sets
function Base.show(io::IO, lc::LocClosedSet)
    string_rep = "V("
    for (i, f) in enumerate(lc.eqns)
        if i != length(lc.eqns)
            string_rep *= "$f, "
        else
            string_rep *= "$(f)) \\ "
        end
    end
    string_rep *= "V("
    for (i, f) in enumerate(lc.ineqns)
        if i != length(lc.ineqns)
            string_rep *= "($f)*"
        else
            string_rep *= "($f)"
        end
    end
    string_rep *= ")"
    print(io, string_rep)
end

function ring(X::LocClosedSet)
    return parent(first(X.eqns))
end

function codim(X::LocClosedSet)
    mis = first(max_ind_sets(X.gb))
    return length(findall(b -> !b, mis))
end

function is_empty_set(X::LocClosedSet)
    R = ring(X)
    if one(R) in X.gb
        return true
    else
        gb2 = saturate(X.gb, last(gens(R)))
        return one(R) in gb2
    end
    return false
end

function add_equation!(X::LocClosedSet, f::MPolyRingElem)
    @info "adding equation"
    push!(X.eqns, f)
    X.gb = saturate(push!(X.gb, f), X.ineqns)
end

function add_inequation!(X::LocClosedSet, h::MPolyRingElem)
    @info "adding inequation"
    push!(X.ineqns, h)
    X.gb = saturate(X.gb, h)
end

function add_inequation(X::LocClosedSet, h::MPolyRingElem)
    if isone(h)
        return X
    end
    Y = deepcopy(X)
    add_inequation!(Y, h)
    return Y
end

function hull(X::LocClosedSet, g::MPolyRingElem)
    gb = X.gb
    sat_gb = saturate(gb, g)
    H = my_normal_form(sat_gb, gb)
    filter!(h -> !iszero(h), H)
    sort!(H, by = h -> total_degree(h))
    return remove(X, H)
end

function remove(X::LocClosedSet,
                H::Vector{<:MPolyRingElem})

    res = typeof(X)[]
    isempty(H) && return res
    h = first(H)
    cells = remove(X, H[2:end])
    if iszero(my_normal_form([h], X.gb))
        return cells
    else
        Y = add_inequation(X, h)
        push!(res, Y)
        for Z in cells
            cells2 = hull(Z, h)
            append!(res, cells2)
        end
        return res
    end
    return res
end

# --- helper functions --- #

function saturate(F::Vector{P}, nzs::Vector{P}) where {P <: MPolyRingElem}
    res = F
    if isempty(nzs)
        res = groebner_basis(Ideal(res), complete_reduction = true)
    else
        for h in nzs
            res = saturate(res, h)
        end
    end
    return res
end

function saturate(F::Vector{P}, nz::P) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), pushfirst!(["x$i" for i in 1:nvars(R)], "t"),
                             ordering = :degrevlex)
    Fconv = typeof(first(F))[]
    for f in F
        ctx = MPolyBuildCtx(S)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            enew = pushfirst!(e, 0)
            push_term!(ctx, c, e)
        end
        push!(Fconv, finish(ctx))
    end

    ctx = MPolyBuildCtx(S)
    for (e, c) in zip(exponent_vectors(nz), coefficients(nz))
        enew = pushfirst!(e, 0)
        push_term!(ctx, c, e)
    end
    push!(Fconv, first(vars)*finish(ctx) - 1)

    gb = groebner_basis(Ideal(Fconv), complete_reduction = true, eliminate = 1)

    # convert back to original ring
    res = Vector{P}(undef, length(gb))

    for (i, p) in enumerate(gb)
        ctx = MPolyBuildCtx(R)
        for (e, c) in zip(exponent_vectors(p), coefficients(p))
            push_term!(ctx, c, e[2:end])
        end
        res[i] = finish(ctx)
    end
    return res
end

function max_ind_sets(gb::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(gb))
    res = [trues(ngens(R))]

    lms = (Nemo.leading_monomial).(gb)
    for lm in lms
        to_del = Int[]
        new_miss = BitVector[]
        for (i, mis) in enumerate(res)
            nz_exps_inds = findall(e -> !iszero(e),
                                   first(Nemo.exponent_vectors(lm)))
            ind_var_inds = findall(mis)
            if issubset(nz_exps_inds, ind_var_inds)
                for j in nz_exps_inds
                    new_mis = copy(mis)
                    new_mis[j] = false
                    push!(new_miss, new_mis)
                end
                push!(to_del, i)
            end
        end
        deleteat!(res, to_del)
        append!(res, new_miss)
        unique!(res)
    end

    max_length = maximum(mis -> length(all(mis)), res)
    filter!(mis -> length(mis) != max_length, res)
    return res
end
    
