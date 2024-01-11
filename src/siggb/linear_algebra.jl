function echelonize!(matrix::MacaulayMatrix,
                     tags::Tags,
                     char::Val{Char},
                     shift::Val{Shift};
                     trace::Bool=true) where {Char, Shift}

    arit_ops = 0

    col2hash = matrix.col2hash
    buffer = zeros(Cbuf, matrix.ncols)
    hash2col = Vector{MonIdx}(undef, matrix.ncols)
    rev_sigorder = Vector{Int}(undef, matrix.nrows)
    pivots = matrix.pivots

    tr_diagonal = if trace
        Vector{Coeff}(undef, matrix.nrows)
    else
        Coeff[]
    end

    tr_mat_data = if trace
        Vector{Vector{Tuple{Int, Coeff}}}(undef, matrix.nrows)
    else
        Vector{Tuple{Int, Coeff}}[]
    end

    tr_mat = TracerMatrix(Dict{Sig, Tuple{Int, Int}}(),
                          Dict{Int, Sig}(),
                          tr_diagonal,
                          tr_mat_data)

    @inbounds for i in 1:matrix.nrows
        rev_sigorder[matrix.sig_order[i]] = i
        row_ind = matrix.sig_order[i]
    end

    @inbounds for i in 1:matrix.ncols
        hash2col[col2hash[i]] = MonIdx(i)
    end

    @inbounds for i in 1:matrix.nrows
        row_ind = matrix.sig_order[i]

        # store tracer data
        row_sig = matrix.sigs[row_ind]
        if trace
            tr_mat.rows[row_sig] = (row_ind, matrix.parent_inds[row_ind])

            # allocate a row for the tracer matrix
            # at most we subtract (i-1) other rows
            row_ops = Vector{Tuple{Int, Coeff}}(undef, i - 1)
            tr_mat.col_inds_and_coeffs[row_ind] = row_ops
            tr_mat.row_ind_to_sig[row_ind] = row_sig
            tr_mat.diagonal[row_ind] = one(Coeff)
        end

        row_cols = matrix.rows[row_ind]
        l_col_idx = hash2col[first(row_cols)]
        if pivots[l_col_idx] == row_ind
            trace && resize!(row_ops, 0)
            continue
        end

        # check if the row can be reduced
        does_red = false
        for (j, m_idx) in enumerate(row_cols)
            colidx = hash2col[m_idx]
            pividx = pivots[colidx]
            if !iszero(pividx) && rev_sigorder[pividx] < i
                row_sig_ind = index(row_sig)
                piv_sig_ind = index(matrix.sigs[pividx])
                row_tag = gettag(tags, row_sig_ind)
                piv_tag = gettag(tags, piv_sig_ind)
                if row_tag == :colins
                    does_red = piv_tag == :colins
                    does_red && break
                elseif row_tag == :col
                    does_red = piv_tag != :col
                    does_red && break
                else
                    does_red = true
                    break
                end
            end
        end
        if !does_red
            pivots[l_col_idx] = row_ind
            trace && resize!(row_ops, 0)
            continue
        end

        # buffer the row
        row_coeffs = matrix.coeffs[row_ind]
        @inbounds for (k, j) in enumerate(row_cols)
            col_idx = hash2col[j]
            buffer[col_idx] = row_coeffs[k]
        end
        
        # do the reduction
        n_row_subs = 0
        @inbounds for j in 1:matrix.ncols
            a = buffer[j] % Char
            iszero(a) && continue
            pividx = pivots[j]
            if iszero(pividx) || rev_sigorder[pividx] >= i
                continue
            end

            n_row_subs += 1
            if trace
                row_ops[n_row_subs] = (pividx, a)
            end

            # subtract a*rows[pivots[j]] from buffer
            pivmons = matrix.rows[pividx]
            pivcoeffs = matrix.coeffs[pividx]

            arit_ops_new = critical_loop!(buffer, j, a, hash2col, pivmons,
                                          pivcoeffs, shift)
            arit_ops += arit_ops_new
        end

        # finalize tracer row, add it to tracer matrix
        trace && resize!(row_ops, n_row_subs)

        new_row_length = 0
        @inbounds for j in 1:matrix.ncols
            iszero(buffer[j]) && continue
            buffer[j] = buffer[j] % Char
            iszero(buffer[j]) && continue
            new_row_length += 1
        end

        # write out matrix row again
        j = 1
        inver = one(Coeff)
        new_row = Vector{MonIdx}(undef, new_row_length)
        new_coeffs = Vector{Coeff}(undef, new_row_length)
        @inbounds for k in 1:matrix.ncols
            iszero(buffer[k]) && continue
            if isone(j)
                pivots[k] = row_ind
                inver = inv(Coeff(buffer[k]), char)
            end
            new_row[j] = col2hash[k]
            new_coeffs[j] = isone(j) ? one(Coeff) : mul(inver, buffer[k], char)
            buffer[k] = zero(Cbuf)
            j += 1
        end
        # store that we normalized the row
        if trace
            tr_mat.diagonal[row_ind] = inver
        end

        # check if row lead reduced
        m = monomial(row_sig)
        @inbounds if isempty(new_row) || (matrix.rows[row_ind][1] != new_row[1] && any(!iszero, m.exps))
            # TODO: not super happy with this check
            if !(row_ind in matrix.toadd[1:matrix.toadd_length])
                matrix.toadd[matrix.toadd_length+1] = row_ind
                matrix.toadd_length += 1
            end
        end

        matrix.rows[row_ind] = new_row
        matrix.coeffs[row_ind] = new_coeffs
    end
    if !iszero(arit_ops)
        @info "$(arit_ops) submul's"
    end

    if !is_triangular(matrix)
        print(matrix)
        error("not triangular")
    end

    return tr_mat
end

# subtract mult
# TODO: for module tracking we won't be able to assume that mult = buffer[bufind]
@inline function critical_loop!(buffer::Vector{Cbuf},
                                bufind::Int,
                                mult::Cbuf,
                                hash2col::Vector{MonIdx},
                                pivmons::Vector{MonIdx},
                                pivcoeffs::Vector{Coeff},
                                shift::Val{Shift}) where Shift
    
    @inbounds buffer[bufind] = zero(Cbuf)
    l = length(pivmons)
    @turbo warn_check_args=false for k in 2:l
        c = pivcoeffs[k]
        m_idx = pivmons[k]
        colidx = hash2col[m_idx]
        buffer[colidx] = submul(buffer[colidx], mult, c, shift)
    end
    return l-1
end

# helper functions
# field arithmetic
function maxshift(::Val{Char}) where Char
    bufchar = Cbuf(Char)
    return bufchar << leading_zeros(bufchar)
end

# compute a representation of a - b*c mod char (char ~ Shift)
@inline function submul(a::Cbuf, b::Cbuf, c::Coeff, ::Val{Shift}) where Shift
    r0 = a - b*Cbuf(c)
    r1 = r0 + Shift
    r0 > a ? r1 : r0
end

@inline function inv(a::Coeff, ::Val{Char}) where Char
    return invmod(Cbuf(a), Cbuf(Char)) % Coeff
end

@inline function addinv(a::Coeff, ::Val{Char}) where Char
    return Char - a
end

@inline function mul(a, b, ::Val{Char}) where Char 
    isone(a) && return b
    isone(b) && return a
    return Coeff((Cbuf(a) * Cbuf(b)) % Char)
end

@inline function add(a, b, ::Val{Char}) where Char
    c0 = a + b
    return Coeff(c0 % Char)
    # c1 = c0 - Coeff(Char)
    # return max(c0, c1)
end

# for tracer

function new_tracer()
    mats = TracerMatrix[]
    basis_ind_to_mat = Vector{Int}(undef, init_basis_size)
    syz_ind_to_mat = Vector{Int}(undef, init_syz_size)
    return Tracer(mats, basis_ind_to_mat,
                  syz_ind_to_mat, 0,
                  init_basis_size)
end

# for debug helping

function is_triangular(matrix::MacaulayMatrix)
    lms = [first(row) for row in matrix.rows[1:matrix.nrows] if !isempty(row)]
    return length(lms) == length(unique(lms))
end

function assert_sigs(matrix::MacaulayMatrix)
    @inbounds sgs = matrix.sigs[1:matrix.nrows]
    return length(sgs) == length(unique(sgs))
end

function Base.show(io::IO, matrix::MacaulayMatrix)
    print_mat = zeros(Int, matrix.nrows, matrix.ncols)
    col2hash = matrix.col2hash
    hash2col = Vector{MonIdx}(undef, matrix.ncols)
    @inbounds for i in 1:matrix.ncols
        hash2col[col2hash[i]] = MonIdx(i)
    end
    for i in 1:matrix.nrows
        row_ind = matrix.sig_order[i]
        # buffer the row
        row_coeffs = matrix.coeffs[row_ind]
        row_cols = matrix.rows[row_ind]
        @inbounds for (k, j) in enumerate(row_cols)
            col_idx = hash2col[j]
            print_mat[i, col_idx] = row_coeffs[k]
        end
    end
    display(print_mat)
end
