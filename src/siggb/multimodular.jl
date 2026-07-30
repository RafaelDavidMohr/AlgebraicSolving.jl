lift_to_int(a::FqFieldElem) = Int(lift(ZZ, a))

# TODO:adjust
function does_not_match(pr::ReconstructPol, p::FqMPolyRingElem)
    ev = collect(exponent_vectors(p))
    return length(ev) != length(pr.exps) || first(ev) != first(pr.exps)
end

function ReconstructPol(p::FqMPolyRingElem)
    exps = collect(exponent_vectors(p))
    mod_coeffs = lift_to_int.(collect(coefficients(p)))
    coeff_cands = (c -> QQ(c)).(mod_coeffs)
    return ReconstructPol(exps, mod_coeffs, coeff_cands, false)
end

# update polynomial at currend index of registry
function update!(reg::ReconstructRegistry,
                 new_pol::FqMPolyRingElem)

    if length(reg.pols) < reg.curr_ind
        push!(reg.pols, ReconstructPol(new_pol))
        reg.curr_ind += 1
        return
    end

    pr = reg.pols[reg.curr_ind]
    does_not_match(pr, new_pol) && error("Bad prime during multi-modular computation.")

    pr.is_stable && return

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
end

    
