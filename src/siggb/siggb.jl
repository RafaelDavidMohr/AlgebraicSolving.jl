# sizes for initialization
const init_ht_size = 2^17
const init_basis_size = 10000
const init_syz_size = 1000
const init_pair_size = 10000
# default sorting alg
const def_sort_alg = Base.Sort.DEFAULT_UNSTABLE
include("typedefs.jl")
include("monomials.jl")
include("hashtable.jl")
include("rewriting.jl")
include("update.jl")
include("symbolic_pp.jl")
include("linear_algebra.jl")
include("module.jl")


#---------------- user functions --------------------#

@doc Markdown.doc"""
    sig_groebner_basis(sys::Vector{T}; info_level::Int = 0, degbound::Int = 0) where {T <: MPolyRingElem}

Compute a Signature Gröbner basis of the sequence `sys` w.r.t. to the
degree reverse lexicographical monomial ordering and the degree
position-over-term ordering induced by `sys`. The output is a vector
of `Tuple{Tuple{Int64, T}, T}` where the first element indicates the
signature and the second the underlying polynomial.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.
**Note**: The input generators must be homogeneous.
**Note**: The algorithms behaviour may depend heavily on how the elements in `sys` are sorted.

# Arguments
- `sys::Vector{T} where T <: MpolyElem`: input generators.
- `info_level::Int=0`: info level printout: off (`0`, default), computational details (`1`)
- `degbound::Int=0`: degree bound for Gröbner basis computation, compute a full Gröbner basis if `0` (default) or only up to degree `d`.

# Example
```jldoctest
julia> using AlgebraicSolving

julia> R, vars = polynomial_ring(GF(17), ["x$i" for i in 1:4])
(Multivariate polynomial ring in 4 variables over GF(17), fpMPolyRingElem[x1, x2, x3, x4])

julia> F = AlgebraicSolving.cyclic(R)
fpMPolyRingElem[x1 + x2 + x3 + x4, x1*x2 + x1*x4 + x2*x3 + x3*x4, x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4, x1*x2*x3*x4 + 16]

julia> Fhom = AlgebraicSolving._homogenize(F.gens)
4-element Vector{fpMPolyRingElem}:
 x1 + x2 + x3 + x4
 x1*x2 + x2*x3 + x1*x4 + x3*x4
 x1*x2*x3 + x1*x2*x4 + x1*x3*x4 + x2*x3*x4
 x1*x2*x3*x4 + 16*x5^4

julia> sig_groebner_basis(Fhom)
7-element Vector{Tuple{Tuple{Int64, fpMPolyRingElem}, fpMPolyRingElem}}:
 ((1, 1), x1 + x2 + x3 + x4)
 ((2, 1), x2^2 + 2*x2*x4 + x4^2)
 ((3, 1), x2*x3^2 + x3^2*x4 + 16*x2*x4^2 + 16*x4^3)
 ((4, 1), x2*x3*x4^2 + x3^2*x4^2 + 16*x2*x4^3 + x3*x4^3 + 16*x4^4 + 16*x5^4)
 ((4, x3), x3^3*x4^2 + x3^2*x4^3 + 16*x3*x5^4 + 16*x4*x5^4)
 ((4, x2), x2*x4^4 + x4^5 + 16*x2*x5^4 + 16*x4*x5^4)
 ((4, x2*x3), x3^2*x4^4 + x2*x3*x5^4 + 16*x2*x4*x5^4 + x3*x4*x5^4 + 15*x4^2*x5^4)
```
"""
function sig_groebner_basis(sys::Vector{T}; info_level::Int=0, degbound::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)
    
    sysl = length(sys)

    # fill basis and pairset
    basis, pairset, tags = fill_data_structs(sys_mons, sys_coeffs, basis_ht, sysl+1)

    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    with_logger(logger) do
        siggb!(basis, pairset, basis_ht, char, shift, tags,
               degbound = degbound)
    end

    # output
    R = parent(first(sys))
    eltp = typeof(first(sys))
    outp = Tuple{Tuple{Int, eltp}, eltp}[]
    @inbounds for i in basis.basis_offset:basis.basis_load
        pol = convert_to_pol(R,
                             [basis_ht.exponents[m] for m in basis.monomials[i]],
                             basis.coefficients[i])

        s = basis.sigs[i]
        ctx = MPolyBuildCtx(R)
        push_term!(ctx, 1, Vector{Int}(monomial(s).exps))
        sig = (Int(index(s)), finish(ctx))

        push!(outp, (sig, pol))
    end

    return outp
end

function sig_decomp(sys::Vector{T}; info_level::Int=0) where {T <: MPolyRingElem}

    # data structure setup/conversion
    sys_mons, sys_coeffs, basis_ht, char, shift = input_setup(sys)
    
    sysl = length(sys)

    # fill basis and pairset
    basis, pairset, tags = fill_data_structs(sys_mons, sys_coeffs,
                                             basis_ht, sysl+1, :split)

    # compute divmasks
    fill_divmask!(basis_ht)
    @inbounds for i in 1:sysl
        basis.lm_masks[i] = basis_ht.hashdata[basis.monomials[i][1]].divmask
    end

    logger = ConsoleLogger(stdout, info_level == 0 ? Warn : Info)
    result = with_logger(logger) do
        sig_decomp!(basis, pairset, basis_ht, char, shift, tags)
    end

    # output
    R = parent(first(sys))
    eltp = typeof(first(sys))
    outp = LocClosedSet{eltp}[]
    for i in 1:length(result)
        lc_set = LocClosedSet(eltp[], eltp[])
        bs, tgs = result[i]
        @inbounds for j in 1:bs.input_load
            bs.is_red[j] && continue
            s_ind = index(bs.sigs[j])
            pol = convert_to_pol(R,
                                 [basis_ht.exponents[m] for m in bs.monomials[j]],
                                 bs.coefficients[j])
            if gettag(tgs, s_ind) == :split
                push!(lc_set.eqns, pol)
            elseif gettag(tgs, s_ind) == :col
                push!(lc_set.ineqns, pol)
            end
        end
        push!(outp, lc_set)
    end

    return outp
end

#---------------- function for sig_groebner_basis --------------------#

function siggb!(basis::Basis{N},
                pairset::Pairset,
                basis_ht::MonomialHashtable,
                char::Val{Char},
                shift::Val{Shift},
                tags::Tags;
                degbound = 0) where {N, Char, Shift}

    # index order
    ind_order = IndOrder((SigIndex).(collect(1:basis.basis_offset-1)),
                         SigIndex(basis.basis_offset-1))

    # tracer
    tr = new_tracer()

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    while !iszero(pairset.load)
        if !iszero(degbound) && first(pairset.elems).deg > degbound
            break
        end
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        deg = select_normal!(pairset, basis, matrix,
                             basis_ht, symbol_ht, ind_order, tags)
        symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                     ind_order, tags)
        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        tr_mat = echelonize!(matrix, tr, tags, char, shift)

        push!(tr.mats, tr_mat)

        update_siggb!(basis, matrix, pairset, symbol_ht,
                      basis_ht, ind_order, tags,
                      tr, char)
        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end
end


#---------------- functions for splitting --------------------#

function sig_decomp!(basis::Basis{N},
                     pairset::Pairset,
                     basis_ht::MonomialHashtable,
                     char::Val{Char},
                     shift::Val{Shift},
                     tags::Tags) where {N, Char, Shift}

    queue = [(basis, pairset, tags)]
    result = Tuple{Basis{N}, Tags}[]

    while !isempty(queue)
        @info "starting component"
        bs, ps, tgs = popfirst!(queue)
        found_zd, zd_coeffs, zd_mons, zd_ind = siggb_for_split!(bs, ps, tgs,
                                                                basis_ht, char,
                                                                shift)
        if found_zd
            @info "splitting component"
            bs1, ps1, tgs1, bs2, ps2, tgs2 = split!(bs, basis_ht, zd_mons,
                                                    zd_coeffs, zd_ind, tgs)
            push!(queue, (bs2, ps2, tgs2))
            push!(queue, (bs1, ps1, tgs1))
        else
            @info "finished component"
            push!(result, (bs, tgs))
        end
        @info "------------------------------------------"
    end
    return result
end

function siggb_for_split!(basis::Basis{N},
                          pairset::Pairset,
                          tags::Tags,
                          basis_ht::MonomialHashtable,
                          char::Val{Char},
                          shift::Val{Shift}) where {N, Char, Shift}

    # index order
    ind_order = IndOrder((SigIndex).(collect(1:basis.basis_offset-1)),
                         SigIndex(basis.basis_offset-1))

    # tracer
    tr = new_tracer()

    # syzygy queue
    syz_queue = Sig{N}[]

    # lms of nonzero conditions
    @inbounds nz_from = findfirst(sig -> gettag(tags, index(sig)) == :col,
                                  basis.sigs[1:basis.input_load])
    nz_lms = [one_monomial(Monomial{N})]
    @inbounds if !isnothing(nz_from)
        for i in nz_from:basis.input_load
            lm_hash_idx = first(basis.monomials[i])
            lm = basis_ht.exponents[lm_hash_idx]
            push!(nz_lms, lm)
        end
    end
    nz_length = length(nz_lms)
    max_nz_deg = maximum(lm -> lm.deg, nz_lms)

    sort_pairset_by_degree!(pairset, 1, pairset.load-1)

    while !iszero(pairset.load)
	matrix = initialize_matrix(Val(N))
        symbol_ht = initialize_secondary_hash_table(basis_ht)

        deg = select_normal!(pairset, basis, matrix,
                             basis_ht, symbol_ht, ind_order, tags)
        symbolic_pp!(basis, matrix, basis_ht, symbol_ht,
                     ind_order, tags)

        finalize_matrix!(matrix, symbol_ht, ind_order)
        iszero(matrix.nrows) && continue
        tr_mat = echelonize!(matrix, tr, tags, char, shift)
        push!(tr.mats, tr_mat)
        tr.deg_to_mat[deg] = length(tr.mats)

        update_siggb!(basis, matrix, pairset, symbol_ht,
                      basis_ht, ind_order, tags,
                      tr, char, syz_queue)

        # check to see if we can split with one of the syzygies
        @inbounds while !isempty(syz_queue)
            @info "checking known syzygies"
            syz_sig = first(syz_queue)
            syz_mon = monomial(syz_sig)
            if syz_mon.deg + max_nz_deg <= deg
                popfirst!(syz_queue)
                
                # membership check with leading monomials
                syz_mon_tms_nz = [mul(syz_mon, lm) for lm in nz_lms]
                syz_masks_tms_nz = [divmask(mon, basis_ht.divmap, basis_ht.ndivbits)
                                    for mon in syz_mon_tms_nz]
                does_div = false
                for j in basis.basis_offset:basis.basis_load
                    lm = basis_ht.exponents[first(basis.monomials[j])]
                    lm_msk = basis.lm_masks[j]
                    for i in 1:nz_length
                        mon = syz_mon_tms_nz[i]
                        msk = syz_masks_tms_nz[i]
                        if divch(lm, mon, lm_msk, msk)
                            does_div = true
                            break
                        end
                    end
                end
                does_div && break

                # from here on we assume that the membership check passed
                @info "membership check passed, constructing module rep"
                syz_ind = index(syz_sig)
                tr_ind = tr.deg_to_mat[syz_mon.deg + basis.degs[syz_ind]]
                cofac_coeffs, cofac_mons = construct_module(syz_sig, basis, tr_ind, tr, char,
                                                            ind_order.max_ind,
                                                            ind_order, syz_ind)[syz_ind]

                # normalize cofac coefficients
                inver = inv(first(cofac_coeffs), char)
                @inbounds for i in eachindex(cofac_coeffs)
                    if isone(i)
                        cofac_coeffs[i] = one(Coeff)
                        continue
                    end
                    cofac_coeffs[i] = mul(inver, cofac_coeffs[i], char)
                end

                cofac_mons_hashed = [insert_in_hash_table!(basis_ht, mon) for mon in cofac_mons]
                return true, cofac_coeffs, cofac_mons_hashed, syz_ind
            else
                break
            end
        end

        sort_pairset_by_degree!(pairset, 1, pairset.load-1)
    end

    return false, Coeff[], MonIdx[], zero(SigIndex)
end

function split!(basis::Basis,
                basis_ht::MonomialHashtable,
                cofac_mons::Vector{MonIdx},
                cofac_coeffs::Vector{Coeff},
                zd_ind::SigIndex,
                tags::Tags)

    @inbounds begin
        # 1st component
        sys1_mons = copy(basis.monomials[1:basis.input_load])
        sys1_coeffs = copy(basis.coefficients[1:basis.input_load])
        sys1l = length(sys1_mons)
        nz_from1 = findfirst(sig -> gettag(tags, index(sig)) == :col,
                            basis.sigs[1:basis.input_load])
        if isnothing(nz_from1)
            nz_from1 = sys1l + 1 
        end
        
        # find out where to insert zero divisor
        zd_deg = basis_ht.exponents[first(cofac_mons)].deg
        ins_ind = findfirst(d -> d > zd_deg, basis.degs)
        if isnothing(ins_ind)
            ins_ind = sys1l + 1
        end

        # insert zd in system
        insert!(sys1_mons, ins_ind, cofac_mons)
        insert!(sys1_coeffs, ins_ind, cofac_coeffs)

        # build basis/pairset/tags for first new system
        basis1, pairset1, tags1 = fill_data_structs(sys1_mons, sys1_coeffs,
                                                    basis_ht, nz_from1 + 1,
                                                    :split)

        # 2nd component
        sys2_mons = basis.monomials[1:basis.input_load]
        sys2_coeffs = basis.coefficients[1:basis.input_load]
        deleteat!(sys2_mons, zd_ind)
        deleteat!(sys2_coeffs, zd_ind)
        sys2l = length(sys2_mons)
        nz_from2 = nz_from1 - 1

        # append zd as nonzero condition
        push!(sys2_mons, copy(cofac_mons))
        push!(sys2_coeffs, copy(cofac_coeffs))

        # build basis/pairset/tags for second new system
        basis2, pairset2, tags2 = fill_data_structs(sys2_mons, sys2_coeffs,
                                                    basis_ht, nz_from2, :split)
    end

    return basis1, pairset1, tags1, basis2, pairset2, tags2
end    


#---------------- functions for setting up data structures --------------------#

function input_setup(sys::Vector{<:MPolyRingElem})

    if isempty(sys)
        error("Input system is empty.")
    end
    
    R = first(sys).parent
    Rchar = characteristic(R)

    # check if input is ok
    if Rchar > 2^31 || iszero(Rchar)
        error("At the moment we only support finite fields up to prime characteristic < 2^31.")
    end
    sysl = length(sys)
    degs = Vector{Exp}(undef, sysl)
    @inbounds for (i, f) in enumerate(sys)
        deg = total_degree(f)
        if deg > typemax(Exp)
            error("input degrees too large.")
        end
        degs[i] = Exp(deg)
        for m in exponent_vectors(f)
            if sum(m) != deg
                error("input system must be homogeneous.")
            end
        end
    end

    # constants for fast arithmetic
    char = Val(Coeff(Rchar.d))
    shift = Val(maxshift(char))

    # convert to and initialize our data structures
    nv = nvars(R)
    basis_ht = initialize_basis_hash_table(Val(nv))

    sys_mons = Vector{Vector{MonIdx}}(undef, sysl)
    sys_coeffs = Vector{Vector{Coeff}}(undef, sysl)

    # store initial pols in basis and pairset
    @inbounds for i in 1:sysl
        f = sys[i]
        lf = length(f)

        # gather up monomials and coeffs
        exps = collect(exponent_vectors(f))
        cfs = collect(coefficients(f))
        mons = Vector{MonIdx}(undef, lf)
        coeffs = Vector{Coeff}(undef, lf)
        inver = one(Coeff)
        @inbounds for j in 1:lf
            m = monomial(SVector{nv}((Exp).(exps[j])))
            eidx = insert_in_hash_table!(basis_ht, m)
            if isone(j)
                inver = inv(Coeff(cfs[1].data), char)
            end
            cf = isone(j) ? one(Coeff) : mul(inver, Coeff(cfs[j].data), char)
            mons[j] = eidx
            coeffs[j] = cf
        end
        s = sortperm(mons, by = eidx -> basis_ht.exponents[eidx],
                     lt = lt_drl, rev = true)
        mons = mons[s]
        coeffs = coeffs[s]
        sys_mons[i] = copy(mons)
        sys_coeffs[i] = copy(coeffs)
    end

    return sys_mons, sys_coeffs, basis_ht, char, shift
end


function fill_data_structs(sys_mons::Vector{Vector{MonIdx}},
                           sys_coeffs::Vector{Vector{Coeff}},
                           basis_ht::MonomialHashtable{N},
                           nz_from::Int,
                           def_tag::Symbol=:seq) where N

    # initialize basis
    sysl = length(sys_mons)
    basis = new_basis(init_basis_size, init_syz_size, sysl, Val(N))

    # initialize pairset
    pairset = Pairset{N}(Vector{SPair{N}}(undef, init_pair_size),
                         0,
                         init_pair_size)

    # tags
    tags = Tags()

    @inbounds for i in 1:sysl
        s_ind = SigIndex(i)
        if i >= nz_from
            tags[s_ind] = :col
        elseif def_tag != :seq
            tags[s_ind] = :split
        end
        mons = sys_mons[i]
        coeffs = sys_coeffs[i]
        lm = basis_ht.exponents[first(mons)]
        lm_mask = divmask(lm, basis_ht.divmap, basis_ht.ndivbits)

        add_input_element!(basis, pairset, s_ind,
                           sys_mons[i], sys_coeffs[i],
                           lm_mask, lm)
    end

    return basis, pairset, tags
end

function convert_to_pol(R::MPolyRing,
                        exps::Vector{<:Monomial},
                        coeffs::Vector{Coeff})

    ctx = MPolyBuildCtx(R)
    for (e, c) in zip(exps, coeffs)
        push_term!(ctx, c, Vector{Int}(e.exps))
    end
    return finish(ctx)
end

function new_basis(basis_size, syz_size,
                   input_length, ::Val{N}) where N

    sigs = Vector{Sig{N}}(undef, basis_size)
    sigmasks = Vector{MaskSig}(undef, basis_size)
    sigratios = Vector{Monomial{N}}(undef, basis_size)
    rewrite_nodes = Vector{Vector{Int}}(undef, basis_size+1)
    lm_masks = Vector{DivMask}(undef, basis_size)
    monomials = Vector{Vector{MonIdx}}(undef, basis_size)
    coeffs = Vector{Vector{Coeff}}(undef, basis_size)
    is_red = Vector{Bool}(undef, basis_size)
    syz_sigs = Vector{Monomial{N}}(undef, syz_size)
    syz_masks = Vector{MaskSig}(undef, syz_size)
    basis = Basis(sigs, sigmasks, sigratios, rewrite_nodes,
                  lm_masks, monomials, coeffs, is_red,
                  syz_sigs, syz_masks, Exp[],
                  input_length,
                  init_basis_size, 0, input_length,
                  input_length + 1, 0,
                  init_syz_size)

    # root node
    basis.rewrite_nodes[1] = [-1, -1]

    return basis
end

function add_input_element!(basis::Basis{N},
                            pairset::Pairset,
                            ind::SigIndex,
                            mons::Vector{MonIdx},
                            coeffs::Vector{Coeff},
                            lm_divm::DivMask,
                            lm::Monomial) where N

    @inbounds begin
        one_mon = one_monomial(Monomial{N})
        zero_sig = (zero(SigIndex), one_mon)

        # signature
        sig = (ind, one_mon)

        l = basis.input_load + 1

        # store stuff in basis
        basis.sigs[l] = sig
        basis.sigmasks[l] = (ind, zero(DivMask))
        basis.sigratios[l] = lm
        basis.rewrite_nodes[l+1] = [-1, 1]
        basis.monomials[l] = mons
        basis.coefficients[l] = coeffs
        basis.is_red[l] = false
        push!(basis.degs, lm.deg)
        basis.lm_masks[l] = lm_divm
        basis.input_load += 1

        # add child to rewrite root
        push!(basis.rewrite_nodes[1], l+1)
        basis.rewrite_nodes[1][1] += 1

        # add unitvector as pair
        pairset.elems[pairset.load+1] = SPair{N}(sig, zero_sig, zero(DivMask),
                                                 zero(DivMask), Int(ind),
                                                 0, lm.deg)
        pairset.load += 1
    end
end


#---------------- helper functions --------------------#

# write stuff from index i in basis1 to index j in basis2
function overwrite!(basis1::Basis,
                    basis2::Basis,
                    i::Int, j::Int)

    @inbounds begin
        basis2.sigs[j]          = basis1.sigs[i]
        basis2.sigmasks[j]      = basis1.sigmasks[i]
        basis2.sigratios[j]     = basis1.sigratios[i]

        rnodes = copy(basis1.rewrite_nodes[i+1])
        basis2.rewrite_nodes[j+1] = rnodes

        basis2.lm_masks[j]      = basis1.lm_masks[i]

        # TODO: this does not copy, is that safe?
        basis2.monomials[j]     = basis1.monomials[i]
        basis2.coefficients[j]  = basis1.coefficients[i]
        basis2.is_red[j]        = basis1.is_red[i]
    end
end

# homogenize w.r.t. the last variable
function _homogenize(F::Vector{P}) where {P <: MPolyRingElem}
    R = parent(first(F))
    S, vars = polynomial_ring(base_ring(R), ["x$i" for i in 1:nvars(R)+1],
                             ordering = :degrevlex)
    res = typeof(first(F))[]
    for f in F
        ctx = MPolyBuildCtx(S)
        d = total_degree(f)
        for (e, c) in zip(exponent_vectors(f), coefficients(f))
            enew = push!(e, d - sum(e))
            push_term!(ctx, c, e)
        end
        push!(res, finish(ctx))
    end
    return res
end

# test against msolve
function _is_gb(gb::Vector{Tuple{Tuple{Int, P}, P}}) where {P <: MPolyRingElem}
    gb_pols = [p[2] for p in gb]
    gb_msolve = groebner_basis(Ideal(gb_pols), complete_reduction = true)
    
    lms_gb = (Nemo.leading_monomial).(gb_pols)
    lms_msolve = (Nemo.leading_monomial).(gb_msolve)
    res1 = all(u -> any(v -> divides(u, v)[1], lms_gb), lms_msolve)
    res2 = all(u -> any(v -> divides(u, v)[1], lms_msolve), lms_gb)
    return res1 && res2
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
