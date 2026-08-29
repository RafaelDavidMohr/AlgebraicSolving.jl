```@meta
CurrentModule = AlgebraicSolving
DocTestSetup = quote
  using AlgebraicSolving
end
```

```@setup algebraicsolving
using AlgebraicSolving
```

```@contents
Pages = ["real-algebraic-geometry.md"]
```

# Real Algebraic Geometry

## Introduction

AlgebraicSolving allows to work with real semi-algebraic sets, i.e. subsets of real affine space
defined by polynomial equations and inequalities.

## Functionality

AlgebraicSolving implements finding sample points per connected
components of semi-algebraic sets:

```@docs
    points_per_components(
		eqs::Vector{QQMPolyRingElem},
		pos::Vector{QQMPolyRingElem},
		ineqs::Vector{QQMPolyRingElem};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		info_level::Int=0
		)
```

AlgebraicSolving also implements functionality for real root classification:

```@docs
    real_root_classification(
		I::ParametricIdeal{K} where {K<:FracFieldElem},
		g::Vector{MPoly{K}} where {K<:FracFieldElem};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		info_level::Int=0,
		show_progress::Bool=info_level >= 1,
		output_form::Symbol=:sign,
		ignore_no_real_roots::Bool=false,
		)
```

This is based on the computation of Hermite matrices:

```@docs
    hermite_matrix(gb::Vector{T}, b::Vector{T}, g::T) where {T<:MPolyRingElem}
```

```@docs
    hermite_matrix(
		I::ParametricIdeal{K},
		g::MPoly{K};
		retry::Int=10,
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```

```@docs
	hermite_matrix(
		g::Vector{MPoly{K}},
		α::Vector{Int},
		I::ParametricIdeal{K};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```

```@docs
	hermite_matrix(
		g::Vector{MPoly{K}},
		α::Vector{Int},
		I::ParametricIdeal{K},
		vals::Vector{QQFieldElem};
		nr_thrds::Int=1,
		worker_pool::AbstractWorkerPool=default_worker_pool(),
		show_progress::Bool=false
		) where {K<:FracFieldElem}
```

AlgebraicSolving allows to compute a roadmap for the real trace of the zero-set
of the ideal spanned by given input generators over the rationals.

It assumes that the underlying algebraic set is **smooth**, and its real trace is **bounded**.

The underlying engine is provided by msolve.

```@docs
    roadmap(
      I::Ideal{P} where P <: QQMPolyRingElem;
      C::Vector{Vector{Vector{QQFieldElem}}}=Vector{Vector{QQFieldElem}}[],
      info_level::Int=0
    )
```

In addition, AlgebraicSolving can compute equations definition critical loci of polynomial maps over the given algebraic set.

```@docs
   computepolar(
        J::Union{Vector{Int},UnitRange{Int}},
        V::Ideal{P};
        phi::Vector{P} = P[],
        dimproj = length(J)-1,
        only_mins = false
    ) where (P <: MPolyRingElem)
```

### Algorithms for connectivity queries on real algebraic curve

#### Introduction
AlgebraicSolving allows to compute a graph whose associated (3D) piecewise linear curve is semi-algebraically homeomorphic to a real algebraic curve given as input.
This input can be:
  * an `Ideal` given by generators over the rationals, whose real trace is the curve;
  * a one-dimensionnal parametrization given as a `RationalCurveParametrization`. We refer to the `Types` page for such a structure.

It is also possible to compute the arrangement of several curves by simply providing the corresponding vector of Ideals.

In both cases, user can, in addition, provide control points that will be correctly identified in the final graph, from a topological point of view.

The underlying engines are
  * internal subresultants computations for paramtrization of critical points;
  * `msolve` for univariate real root isolation.

#### Computing curve graphs

```@docs
    curve_graph(
        I::Ideal{P}, args...;
        generic=true,
        outf=true,
        v=0,
        kwargs...
      ) where {P <: QQMPolyRingElem}
```

```@docs
    curve_arrangement_graph(
        curves::Vector{Ideal{P}};
        generic=true,
        outf=true,
        v=0,
        kwargs...
    ) where {P <: QQMPolyRingElem}
```

#### Computing with curve graphs

Once such a graph is computed, it is encoded in a CurveGraph object, whose description can be found in the Types section of this documentation.
These are nothing but sets of vertices and nodes.
One can perform different operations on such strutures.

```@docs
    connected_components(
        G::CurveGraph{T}
    ) where T
```

```@docs
    group_by_component(
        G::CurveGraph
    )
```

```@docs
    number_of_connected_components(
        G::CurveGraph
    )
```

```@docs
    merge_graphs(
        graphs::Vector{CurveGraph}
    )
```

#### Displaying curve graph

As e.g. Plots.jl is not part of AlgebraicSolving.jl, it cannot provide direct functions for plotting curve graphs. However, functions for exporting curve graphs in a format adapted to most ploting libraries, encoded in a `GraphPlotData` object. We refer to the Types section for more information on this structuree.

```@docs
    build_graph_data(
        G::CurveGraph{T};
        width=3.0,
        vert=true,
        color="#FF0000"
    )  where T <: Union{Float64,QQFieldElem}
```

```@docs
    build_graph_data(
        CG::Vector{CurveGraph{T}};
        width=3.0,
        vert=true
    )  where T <: Union{Float64,QQFieldElem}
```