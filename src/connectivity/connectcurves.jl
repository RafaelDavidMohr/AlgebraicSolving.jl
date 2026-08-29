include("tools.jl")
include("datastruct.jl")
include("buildpoly.jl")
include("paramcrit.jl")
include("arbtools.jl")
include("isolateboxes.jl")
include("graph.jl")
include("plots.jl")

@doc Markdown.doc"""
    curve_arrangement_graph(curves::Vector{Ideal}; generic::Union{ZZRingElem, Int64}=0, outf=true, v=0, kwargs...)

Computes the combined planar graph of an arrangement of multiple space curves in complete intersection.
Automatically computes the mutual intersections between all curves and guarantees
they are projected using a shared, unified linear form.
"""
function curve_arrangement_graph(curves::Vector{Ideal{P}}, Icon=nothing; generic::Union{ZZRingElem, Int64} = 0, outf=true, v=0, kwargs...) where {P <: QQMPolyRingElem}
    N = length(curves)
    @assert N > 0 "Must provide at least one curve."

    R = parent(curves[1])
    n = nvars(R)
    typeout = outf ? Float64 : QQFieldElem


    # 1. Establish the Shared Projection Context
    # We MUST generate the linear forms here and force all curves/intersections to use them.
    if generic != 0
        Bgen = ZZ(generic < 0 ? 100 : generic)
        make_vec(i) = [j <= n ? rand(-Bgen:Bgen) : (j == n+i ? one(ZZRingElem) : zero(ZZRingElem)) for j in 1:n+3]
        cfs_lfs = make_vec.(1:3)
    else
        cfs_lfs = [ [ j == n-3+i ? -one(ZZRingElem) : (j == n+i ? one(ZZRingElem) : zero(ZZRingElem)) for j in 1:n+3 ] for i in 1:3 ]
    end

    # 3. Compute Intersections and Build Individual Graphs
    graphs = Vector{CurveGraph{typeout}}(undef, N)
    p_inter = Dict{Set{Int}, RationalParametrization}()

    for i in 1:N
        v > 0 && println("Compute graph of curve number $i/$N...")
        p_I = curve_rational_parametrization(curves[i], cfs_lfs=cfs_lfs, check_cfs=false)

        # Control pts are intersections with all other curves; use Dict to avoid re-computation
        C_i = Dict{Int, RationalParametrization}( j => p_inter[Set((i,j))] for j in 1:i-1)
        if p_I.elim == -1
            one_x = one(polynomial_ring(QQ,"x")[1])
            for j in i+1:N
                p_inter[Set((i,j))] =
                    RationalParametrization(Symbol[], ZZRingElem[], -one_x, -one_x, QQPolyRingElem[])
            end
        else
            for j in i+1:N
                I_ij = vcat(curves[i].gens, curves[j].gens) |> Ideal
                try
                    C_i[j] = param_newvars(I_ij, p_I.vars, cfs_lfs)
                catch e
                    if e isa ErrorException &&
                        e.msg == "Dimension of ideal is greater than zero, no solutions provided."
                        error("The curves $i and $j are not in complete intersection")
                    end
                    rethrow()
                end
                p_inter[Set((i,j))] =
                    RationalParametrization(C_i[j].vars, C_i[j].cfs_lf, C_i[j].elim, C_i[j].denom, C_i[j].param)
            end

            # Attach User Query / Control Points with Non-Positive Keys (<= 0)
            if !isnothing(Icon)
                Icons = Icon isa Vector ? Icon : [Icon]
                for (k, Con) in enumerate(Icons)
                    C_i[-k] = param_newvars(Ideal(vcat(curves[i].gens, Con.gens)), p_I.vars, cfs_lfs)
                end
            end
        end

        graphs[i] = curve_graph(p_I, C_i; outf=outf, v=v-1, kwargs...)
    end

    # 4. Merge the graph according to their intersections
    return merge_graphs(graphs)
end

@doc Markdown.doc"""
    curve_graph(I::Ideal; generic=true, outf=true, kwargs...)
    curve_graph(P::CurveRationalParametrization; kwargs...)

    # Both functions accept optional control points C in various formats:
    curve_graph(I, C::Vector{P}; kwargs...)
    curve_graph(I, C::Vector{Vector{P}}; kwargs...)
    curve_graph(I, C::Dict{Int, Vector{P}}; kwargs...)

Computes a planar straight-line graph that is homeomorphic to a real algebraic (space) curve.

### Core Workflow & Pre-processing
1. **Parametrization:** If given an `Ideal`, it first computes a `CurveRationalParametrization`.
    If `generic!=0`, it applies a random linear transformation to place the curve into generic position. The absolute size of the used coefficients is generic if positive, 100 else.
2. **Coefficient Extraction:** Pulls the planar projection `f(x,y) = 0` and the vertical lift `z = g(x,y) / df/dy(x,y)` from the parametrization.
3. **Graph Construction:** Computes bounding boxes for critical points, routes connections,
    and identifies "apparent singularities" (2D crossings that do not intersect in 3D).

### Output Data Structure
Returns a `CurveGraph{T}` object (where `T` is determined by the `outf` flag), containing:
* `Vert::Vector{Tuple{T, T}}`: 2D coordinates of the graph vertices (critical points, routing nodes, control points).
* `Edg::Vector{Tuple{Int, Int}}`: Index pairs defining undirected edges between `Vert` indices.
* `Vcon::Dict{Int, Vector{Int}}`: A mapping from the original index of a control point keys its vertex index in `Vert`.

### Arguments
* **`I`** (`Ideal`): The algebraic ideal defining the curve.
* **`P`** (`CurveRationalParametrization`): A pre-computed rational parametrization.
* **`C`** (Optional): User-defined plane points on the curve (control points). Given either as Ideal or RationalParametrization.
* **`generic`** (`Union{ZZRingElem, Int64}`, default `0`): If `!=0`, applies a random shear transformation.
* **`precx`** (`Int`, default `150`): Base numerical precision for real root isolation.
* **`v`** (`Int`, default `0`): Verbosity level.
* **`force_app`** (`Bool`, default `false`): Skips 3D intersection checks, treating all 2D nodes as apparent singularities.
* **`outf`** (`Bool`, default `true`): Output coordinates as `Float64`. If `false`, outputs exact `QQFieldElem`.

### Example
```jldoctest
julia> using AlgebraicSolving

julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> f = -42*x^2 + 101*x*y - 8*x*z - 88*x + 53*y^2 + 71*y*z + 154*y + 53*z^2 + 2*z - 32;
       g = -7*x^2 - 144*x*y - 126*x*z - 7*x + 30*y^2 - 34*y*z + 72*y + 71*z^2 + 46*z + 5;

julia> I = Ideal([f, g]);

julia> G = curve_graph(I);

julia> number_of_connected_components(G)
3
```
"""
function curve_graph(I::Ideal{P}, args...; generic::Union{ZZRingElem, Int64} = 0, outf=true, v=0, kwargs...) where {P <: QQMPolyRingElem}
    R = parent(I)
    n = nvars(R)
    typeout = outf ? Float64 : QQFieldElem

    # Base Case: 1D Ideal (Points)
    if n == 1
        if dimension(I) == 1
            Vert = [(zero(typeout), zero(typeout))]
        else
            sols = real_solutions(I)
            Vert = [(typeout(xi[1]), zero(typeout)) for xi in sols]
        end
        # If C is a Dict, instantiate keys. If it's a Vector or absent, return empty Dict
        C_keys = length(args) > 0 && args[1] isa Dict ? keys(args[1]) : Int[]
        Vcon = Dict{Int, Vector{Int}}(k => Int[] for k in C_keys)
        return CurveGraph{typeout}(Vert, Tuple{Int, Int}[], Vcon)
    end

    # Pre-processing: Generic Position Shear
    # We explicitly define linear forms if generic != 0 to share them with C
    v > 0 && println("Compute curve rational parametrization...")
    lfs = nothing
    if generic != 0
        Bgen = ZZ(generic < 0 ? ZZ(100) : generic)
        make_vec(i) = [j <= n ? rand(-Bgen:Bgen) : (j == n+i ? one(ZZRingElem) : zero(ZZRingElem)) for j in 1:n+3]
        lfs = make_vec.(1:3)
        u_lfs = lfs[2][1:end-3] # for zero-dim param of control pts
        p_I = curve_rational_parametrization(I, cfs_lfs = lfs, check_cfs=false)
    else
        p_I = curve_rational_parametrization(I)
    end

    # 2. Process Control Ideals (C) if provided
    if length(args) > 0
        C = args[1]
        new_RS = symbols(parent(p_I.elim))

        # Map C based on its input structure
        if C isa AbstractVector
            C_param =  [ isnothing(lfs) ? rational_parametrization(c) : param_use_lfs(c, u_lfs, new_RS[end-1]) for c in C ]
        elseif C isa AbstractDict
            C_param = Dict( k => isnothing(lfs) ? rational_parametrization(c) : param_use_lfs(c, u_lfs, new_RS[end-1]) for (k, c) in C_param)
        else
           error("Control points C must be a Vector or Dict of Ideals.")
        end

        return curve_graph(p_I, C_param; outf=outf, v=v, kwargs...)
    end

    # Delegate to the parametrization function, blindly passing any `C` format down
    return curve_graph(p_I, args...; outf=outf, v=v, kwargs...)
end

# =========================================================================
# MULTIPLE DISPATCH WRAPPERS (INTERNAL)
# =========================================================================

# Small helper
function prepare_param(P,Q)
    R = parent(Q)
    param = is_constant(P.elim) ? [P.elim,P.elim,P.elim] : [P.elim, P.param[end], P.denom]

    return [ evaluate(p, gens(R)[1]) for p in param ]
end

# Case 1: No control points
curve_graph(p::CurveRationalParametrization; kwargs...) =
    _compute_graph_core(p.elim, length(p.param) > 0 ? p.param[end] : zero(parent(p.elim)),
                        Dict{Int, Vector{typeof(p.elim)}}(); kwargs...)

# Case 2: C is a Vector of Parametrizations
curve_graph(p::CurveRationalParametrization, C::Vector{RationalParametrization}; kwargs...) =
    _compute_graph_core(p.elim, length(p.param) > 0 ? p.param[end] : zero(parent(p.elim)),
                        Dict( i => prepare_param(c, p.elim)  for (i,c) in enumerate(C)); kwargs...)

# Case 3: C is a Dictionary of Parametrizations
curve_graph(p::CurveRationalParametrization, C::Dict{Int, RationalParametrization}; kwargs...) =
     _compute_graph_core(p.elim, length(p.param) > 0 ? p.param[end] : zero(parent(p.elim)),
                        Dict( i => prepare_param(c, p.elim) for (i,c) in C); kwargs...)

# =========================================================================
# CORE IMPLEMENTATION
# =========================================================================

function _compute_graph_core(f::P, g::P, C::Dict{Int, Vector{P}};
                             precx=150, v=0, force_app=false, outf=true) where {P <: MPolyRingElem}

    @assert !iszero(f) "Input does not define a curve"

    R = parent(f)
    x, y = gens(R)
    typeout = outf ? Float64 : QQFieldElem

    # Empty set
    total_degree(f) == 0 &&
    return CurveGraph(
    Tuple{typeout, typeout}[],
    Tuple{Int, Int}[],
    Dict{Int, Vector{Int}}()
    )

    # Pre-processing the input
    f, g = int_coeffs([f, g])
    precx = max(2, precx)
    v > 2 && println("f = $f \n g = $g\n $C")

    # Zero-dim param conditions
    d = total_degree(f)
    @assert degree(f, 2) == d && degree(f, 1) == d && degree(g, 2) < d &&
    is_squarefree(f) "Curve not in generic position. Try with generic=false."

    v > 0 && println("Compute parametrization of critical pts...")
    @iftime (v > 0) params = param_crit_split(f, g, v=v-1, force_app=force_app)
    keys_C = collect(keys(C))
    for (i, k) in enumerate(keys_C)
        params[-i] = [ [int_coeffs(C[k][1])], C[k][2], C[k][3] ]
    end

    v > 0 && println("Computing insulating critical boxes")
    @iftime (v > 0) LBcrit, Lprecx = insulate_crit_boxes(f, params, precx, v=v-1)

    v > 0 && println("Compute intersections with critical boxes..")
    @iftime (v > 0) LPCside, LnPCside = intersect_vertical_boxes(f, params, LBcrit, Lprecx, v=v-1)

    # Critical values and their order
    xcrit = Dict(i => [LBcrit[i][j][1] for j in eachindex(LBcrit[i])] for i in keys(LBcrit))
    xcritpermut = order_permut2d(xcrit)

    # Graph data
    Vert = Tuple{QQFieldElem, QQFieldElem}[] # List of points (x,y)
    Edg = Tuple{Int, Int}[] # List of tuples (idx, idy)
    Vcon = Dict{Int, Vector{Int}}(k => Int[] for k in keys(C)) # Index of control vertices

    # Empty or unbounded real curves without control points
    if isempty(xcrit) ||
        (length(xcrit) == 1 && haskey(xcrit, -1) && isempty(xcrit[-1]))
        Vert = reduce(vcat, [ (typeout(e), typeout(yy[1])) for e=0:1 for yy in isolate_eval(f,1,1-2*e)])
        nV = length(Vert) / 2
        Edg = [ (i, i + nV) for i=1:nV ]
        return CurveGraph{typeout}(Vert, Edg, Vcon)
    end

    Lapp_isolated = Tuple{Int, Int}[]
    Lapp_nodes = Tuple{Int, Int}[]

    # Keep track of connection between all computed points when increasing x
    Corr = Dict(m => [CriticalNodeConnections() for _ in xcrit[m]] for m in keys(xcrit))

    for (ind, (i, j)) in enumerate(xcritpermut)
        i1, j1 = ind > 1 ? xcritpermut[ind - 1] : (0, 0) # previous
        i2, j2 = ind < length(xcritpermut) ? xcritpermut[ind + 1] : (0, 0) # next

        # Current cylinder above xcrit[i][j]
        PCside = LPCside[i][j]
        I2L_points = ind < length(xcritpermut) ? LPCside[i2][j2].left.points : nothing

        # Midpoints of the critical box itself
        xcmid = sum(LBcrit[i][j][1]) / 2
        ycmid = sum(LBcrit[i][j][2]) / 2

        nI_left = PCside.left.indices_inside
        nI_right = PCside.right.indices_inside
        ymincrit = minimum(vcat(nI_left, nI_right, [length(PCside.left.points) + 1]))

        # 1. Left side vertices
        # one-to-one connection with previous right side
        if ind > 1
            for k in 1:length(PCside.left.points)
                push!(Corr[i][j].left, Corr[i1][j1].right[k])
            end
        else
            for k in 1:length(PCside.left.points)
                xval = xcrit[i][j][1]
                yval = sum(PCside.left.points[k]) / 2
                push!(Vert, (yval, xval)) # coordinates are reversed in curve_param
                push!(Corr[i][j].left, length(Vert))
            end
        end

        # 2. Right side vertices
        if ind < length(xcritpermut)
            for k in 1:length(PCside.right.points)
                xval = (xcrit[i][j][2] + xcrit[i2][j2][1]) / 2
                yval = sum(PCside.right.points[k] + I2L_points[k]) / 4
                push!(Vert, (yval, xval))
                push!(Corr[i][j].right, length(Vert))
            end
        else
            for k in 1:length(PCside.right.points)
                xval = xcrit[i][j][2]
                yval = sum(PCside.right.points[k]) / 2
                push!(Vert, (yval, xval))
                push!(Corr[i][j].right, length(Vert))
            end
        end

        # 3. Below the critical point
        for k in 1:ymincrit-1
            yval = sum(PCside.left.points[k] + PCside.right.points[k]) / 4
            push!(Vert, (yval, xcmid))
            push!(Corr[i][j].bottom, length(Vert))
            push!(Edg, (Corr[i][j].left[k], length(Vert)))
            push!(Edg, (length(Vert), Corr[i][j].right[k]))
        end

        # 4. The critical point
        if i == 0 && isempty(nI_left) && isempty(nI_right)
            push!(Lapp_isolated, (i, j))
        elseif i == 0
            push!(Edg, (Corr[i][j].left[nI_left[1]], Corr[i][j].right[nI_right[2]]))
            push!(Edg, (Corr[i][j].left[nI_left[2]], Corr[i][j].right[nI_right[1]]))
            push!(Lapp_nodes, (i, j))
        else
            push!(Vert, (ycmid, xcmid))
            push!(Corr[i][j].center, length(Vert))

            for k in nI_left; push!(Edg, (Corr[i][j].left[k], length(Vert))); end
            for k in nI_right; push!(Edg, (length(Vert), Corr[i][j].right[k])); end

            if i < 0 # Control point
                push!(Vcon[keys_C[-i]], length(Vert))
            end
        end

        # 5. Above the critical point
        for k = (length(PCside.left.points) - length(nI_left) - ymincrit + 1):-1:1
            yval = sum(PCside.left.points[end-k+1] + PCside.right.points[end-k+1]) / 4
            push!(Vert, (yval, xcmid))
            push!(Corr[i][j].top, length(Vert))
            push!(Edg, (Corr[i][j].left[end-k+1], length(Vert)))
            push!(Edg, (length(Vert), Corr[i][j].right[end-k+1]))
        end
    end



    if outf
        Vert = [ (Float64(x), Float64(y)) for (x,y) in Vert ]
    end

    # Return the unified data structure
    return CurveGraph{typeout}(Vert, Edg, Vcon)
end