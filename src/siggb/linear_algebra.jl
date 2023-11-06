function echelonize!(matrix::MacaulayMatrix,
                     char::Val{Char},
                     shift::Val{Shift}) where {Char, Shift}

    arit_ops = 0

    col2hash = matrix.col2hash
    buffer = zeros(Cbuf, matrix.ncols)
    cofac_buffer = zeros(Cbuf, matrix.ncofaccols)
    hash2col = Vector{MonIdx}(undef, matrix.ncols)
    rev_sigorder = Vector{Int}(undef, matrix.nrows)
    pivots = matrix.pivots

    @inbounds for i in 1:matrix.nrows
        rev_sigorder[matrix.sig_order[i]] = i
        row_ind = matrix.sig_order[i]
    end

    @inbounds for i in 1:matrix.ncols
        hash2col[col2hash[i]] = MonIdx(i)
    end

    @inbounds for i in 1:matrix.nrows
        row_ind = matrix.sig_order[i]

        row_cols = matrix.rows[row_ind]
        l_col_idx = hash2col[first(row_cols)]
        pivots[l_col_idx] == row_ind && continue

        # check if the row can be reduced
        does_red = false
        for (j, m_idx) in enumerate(row_cols)
            colidx = hash2col[m_idx]
            pividx = pivots[colidx]
            does_red = !iszero(pividx) && rev_sigorder[pividx] < i
            does_red && break
        end
        if !does_red
            pivots[l_col_idx] = row_ind
            continue
        end

        # buffer the row
        row_coeffs = matrix.coeffs[row_ind]
        @inbounds for (k, j) in enumerate(row_cols)
            col_idx = hash2col[j]
            buffer[col_idx] = row_coeffs[k]
        end

        # buffer the cofactor row if applicable
        cofac_ind = get(matrix.row_to_cofac_row, row_ind, zero(row_ind))
        if !iszero(cofac_ind)
            cofac_cols = matrix.cofac_rows[cofac_ind]
            cofac_coeffs = matrix.cofac_coeffs[cofac_ind]
            @inbounds for (k, j) in enumerate(cofac_cols)
                cofac_buffer[j] = cofac_coeffs[k]
            end
        end

        row_sig_ind = index(matrix.sigs[row_ind])
        # do the reduction
        @inbounds for j in 1:matrix.ncols
            a = buffer[j] % Char
            iszero(a) && continue
            pividx = pivots[j]
            if iszero(pividx) || rev_sigorder[pividx] >= i
                continue
            end

            # subtract m*rows[pivots[j]] from buffer
            pivmons = matrix.rows[pividx]
            pivcoeffs = matrix.coeffs[pividx]

            arit_ops_new = critical_loop!(buffer, j, a, hash2col, pivmons,
                                          pivcoeffs, shift)
            arit_ops += arit_ops_new

            # subtract cofactor row if applicable
            piv_cofac_ind = get(matrix.row_to_cofac_row, pividx, zero(row_ind))
            if !iszero(cofac_ind) && !iszero(piv_cofac_ind)
                piv_sig_ind = index(matrix.sigs[pividx])
                if piv_sig_ind == row_sig_ind
                    piv_cofac_mons = matrix.cofac_rows[piv_cofac_ind]
                    piv_cofac_coeffs = matrix.cofac_coeffs[piv_cofac_ind]
                    arit_ops_new = cofac_critical_loop!(cofac_buffer, a, piv_cofac_mons,
                                                        piv_cofac_coeffs, shift)
                    arit_ops += arit_ops_new
                end
            end
        end

        new_row_length = normalize_buffer!(buffer, char, matrix.ncols)

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

        # check if row lead reduced, TODO: dont know if this is reliable
        s = matrix.sigs[row_ind]
        m = monomial(s)
        @inbounds if isempty(new_row) || (matrix.rows[row_ind][1] != new_row[1] && any(!iszero, m.exps))
            matrix.toadd[matrix.toadd_length+1] = row_ind
            matrix.toadd_length += 1
        end

        matrix.rows[row_ind] = new_row
        matrix.coeffs[row_ind] = new_coeffs

        # write out cofactor row if applicable
        if !iszero(cofac_ind)
            j = 1
            new_cofac_row_length = normalize_buffer!(cofac_buffer, char, matrix.ncofaccols)
            new_cofac_row = Vector{MonIdx}(undef, new_cofac_row_length)
            new_cofac_coeffs = Vector{Coeff}(undef, new_cofac_row_length)
            @inbounds for k in 1:matrix.ncofaccols
                iszero(cofac_buffer[k]) && continue
                new_cofac_row[j] = k
                new_cofac_coeffs[j] = mul(inver, cofac_buffer[k], char)
                j += 1
                cofac_buffer[k] = zero(Cbuf)
            end
            matrix.cofac_rows[cofac_ind] = new_cofac_row
            matrix.cofac_coeffs[cofac_ind] = new_cofac_coeffs
        end
    end
    if !iszero(arit_ops)
        @info "$(arit_ops) submul's"
    end
end

# subtract mult
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

@inline function cofac_critical_loop!(buffer::Vector{Cbuf},
                                      mult::Cbuf,
                                      pivmons::Vector{MonIdx},
                                      pivcoeffs::Vector{Coeff},
                                      shift::Val{Shift}) where Shift

    l = length(pivmons)
    @turbo warn_check_args=false for k in 1:l
        c = pivcoeffs[k]
        m_idx = pivmons[k]
        buffer[m_idx] = submul(buffer[m_idx], mult, c, shift)
    end
    return l
end

function normalize_buffer!(buffer::Vector{Cbuf},
                           char::Val{Char},
                           ncols::Int) where Char

    new_row_length = 0
    @inbounds for j in 1:ncols
        iszero(buffer[j]) && continue
        buffer[j] = buffer[j] % Char
        iszero(buffer[j]) && continue
        new_row_length += 1
    end
    return new_row_length
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

@inline function mul(a, b, ::Val{Char}) where Char 
    return Coeff((Cbuf(a) * Cbuf(b)) % Char)
end

# for debug helping

function is_triangular(matrix::MacaulayMatrix)
    lms = [first(row) for row in matrix.rows[1:matrix.nrows] if !isempty(row)]
    return length(lms) == length(unique(lms))
end
