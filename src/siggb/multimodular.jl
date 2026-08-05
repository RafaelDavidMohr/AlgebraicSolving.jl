function reduce_mod_p(f::QQMPolyRingElem, R::FqMPolyRingElem)

    F = base_ring(R)
    ctx = MPolyBuildCtx(R)
    for (cf, exp) in zip(coefficients(f), exponent_vectors(f))
        cfP = F(numerator(cf)) * F(denominator(cf))^(-1)
        push_term!(ctx, cfP, exp)
    end
    return finish(ctx)
end

lift_to_int(a::FqFieldElem) = Int(lift(ZZ, a))

function is_finished(r::ReconstructRegistry)
    !isempty(r.pols) && all(p -> p.is_stable, r.pols)
end

function get_pol(r::ReconstructRegistry, i::Int)
    rp = r.pols[i]
    ctx = MPolyBuildCtx(r.R)
    for (cf, exp) in zip(rp.coeff_cands, rp.exps)
        push_term!(ctx, cf, exp)
    end
    return finish(ctx)
end

function get_pols(r::ReconstructRegistry, inds::Vector{Int})
    return [get_pol(r, i) for i in inds]
end

function is_finished(r::ModularRegistry)
    !isempty(r.pols)
end

function get_pols(r::ModularRegistry, inds::Vector{Int})
    return r.pols[inds]
end

# TODO:adjust
function does_not_match(pr::ReconstructPol, p::FqMPolyRingElem)
    ev = collect(exponent_vectors(p))
    return length(ev) != length(pr.exps) || first(ev) != first(pr.exps)
end

function ReconstructPol(p::FqMPolyRingElem)
    p *= leading_coefficient(p)^(-1)
    exps = collect(exponent_vectors(p))
    mod_coeffs = lift_to_int.(collect(coefficients(p)))
    coeff_cands = (c -> QQ(c)).(mod_coeffs)
    return ReconstructPol(exps, mod_coeffs, coeff_cands, false)
end

# update polynomial at currend index of registry
function update!(reg::ReconstructRegistry,
                 new_pol::FqMPolyRingElem)

    ri = reg.curr_ind
    
    if length(reg.pols) < ri
        push!(reg.pols, ReconstructPol(new_pol))
        reg.curr_ind += 1
        return ri
    end

    pr = reg.pols[ri]
    does_not_match(pr, new_pol) && error("Bad prime during multi-modular computation.")

    if pr.is_stable
        reg.curr_ind += 1
        return ri
    end

    pprod = prod(reg.primes)
    curr_p = reg.current_prime
    all_is_stable = true
    i = 1
    for (ccurr, cnew_fq) in zip(pr.coeff_cands, coefficients(new_pol))
        cnew = lift_to_int(cnew_fq)
        ccurr_new = crt(ZZ(ccurr), ZZ(pprod), ZZ(cnew), ZZ(curr_p))
        pr.mod_coeffs[i] = ccurr_new

        new_qq_coeff = reconstruct(ccurr_new, pprod * curr_p)
        if new_qq_coeff == pr.coeff_cands[i]
            i += 1
            continue
        end
        all_is_stable = false
        pr.coeff_cands[i] = new_qq_coeff

        i += 1
    end
    pr.is_stable = all_is_stable
    reg.curr_ind += 1
    return ri
end

function update!(reg::ModularRegistry{T}, p::T) where T
    push!(reg.pols, p)
    return length(reg.pols)
end
