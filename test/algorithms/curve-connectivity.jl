@testset "Algorithms -> Curve Graph" begin
    A, (x,y)  = polynomial_ring(QQ, [:x,:y])

     # param2 (generic)
    f = -49303382*x^4 + 9395599*x^3*y + 67366686*x^3 - 27407214*x^2*y^2 - 71298014*x^2*y - 24150320*x^2 - 12067817*x*y^3 + 6797370*x*y^2 + 20838420*x*y - 18118613*x + 2357706*y^4 - 13736522*y^3 - 37604516*y^2 - 13868221*y + 6980802
    df = 9395599*x^3 - 54814428*x^2*y - 71298014*x^2 - 36203451*x*y^2 + 13594740*x*y + 20838420*x + 9430824*y^3 - 41209566*y^2 - 75209032*y - 13868221
    g = 32091920*x^4 - 97772598*x^3*y - 245256584*x^3 + 99093410*x^2*y^2 + 162335273*x^2*y + 239310556*x^2 + 28995368*x*y^3 - 59544597*x*y^2 - 20499914*x*y + 98243402*x + 22738226*y^3 + 111105506*y^2 + 86287748*y - 62439535
    P = AlgebraicSolving.CurveRationalParametrization([:z,:x,:y], Vector{Vector{ZZRingElem}}(), f, df, [g])
    G = curve_graph(P)
    @test number_of_connected_components(G) == 3

    # Parametrization of a circle with z=0 but with two new generic variables
    # bug in arb_to_rat fixed
    f = 6829*x^2 - 1482*x*y + 9565*y^2 - 64770304
    df = -1482*x + 19130*y
    g = -1110*x^2 + 11526*x*y + 2076384
    P = AlgebraicSolving.CurveRationalParametrization([:z,:x,:y], Vector{Vector{ZZRingElem}}(), f, df, [g])
    G = curve_graph(P)
    @test number_of_connected_components(G) == 1

    # Tests with ideals
    R,(x1,x2,x3) = polynomial_ring(QQ,[:x1,:x2,:x3])

    # Non-intersecting circle whose projections on (x1,x2) intersect
    I = AlgebraicSolving.Ideal([x3^2 - x3, 2*x1^2 - 2*x1*x3 + 2*x2^2 - 2*x2*x3 + x3 - 2])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 2

    # Same circle, but this time intersecting
    I = AlgebraicSolving.Ideal([x3, 2*x1^4 - 2*x1^3 + 4*x1^2*x2^2 - 2*x1^2*x2 - 3*x1^2 - 2*x1*x2^2 + 2*x1 + 2*x2^4 - 2*x2^3 - 3*x2^2 + 2*x2 + 1])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 1

    ## Plane curves with different topology and especially singularities ##

    # Not generic
    I = AlgebraicSolving.Ideal([x2^2 - x1^2*(x1 + 1), x3])
    @test_throws AssertionError curve_graph(I)
    # Apply change of variables
    I = AlgebraicSolving.Ideal([-x1^3 - 9*x1^2*x2 - x1^2 - 27*x1*x2^2 - 6*x1*x2 - 27*x2^3 - 8*x2^2, x3])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 1

    # Crunode (after generic change of variable)
    f = -2744*x1^3 - 34692*x1^2*x2 + 8085*x1^2 - 146202*x1*x2^2 + 896*x1*x2 - 205379*x2^3 - 3285*x2^2
    I = AlgebraicSolving.Ideal([f, x3])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 1

    # Acnode
    f = 157464*x1^3 + 297432*x1^2*x2 + 7012*x1^2 + 187272*x1*x2^2 + 1496*x1*x2 + 39304*x2^3 + 1445*x2^2
    I = AlgebraicSolving.Ideal([f, x3])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 2

    # Ordinary cusp
    f = 79507*x1^3 + 543606*x1^2*x2 + 9604*x1^2 + 1238916*x1*x2^2 - 3528*x1*x2 + 941192*x2^3 + 324*x2^2
    I = AlgebraicSolving.Ideal([f, x3])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 1

    # Rhamphoid cusp
    f = -243*x1^5 - 25515*x1^4*x2 - 1071630*x1^3*x2^2 - 22504230*x1^2*x2^3 + 8464*x1^2 - 236294415*x1*x2^4 + 18216*x1*x2 - 992436543*x2^5 + 9801*x2^2
    I = AlgebraicSolving.Ideal([f, x3])
    G = curve_graph(I)

    # Non-generic because two crit-point above (but rat_cur_param works)
    f = x1^4 - x1^3 + 2*x1^2*x2^2 - 7//4*x1^2 - x1*x2^2 + x1 + x2^4 - 7//4*x2^2 + 3//4
    I = AlgebraicSolving.Ideal([f, x3])
    @test_throws ErrorException curve_graph(I)
    I = AlgebraicSolving.Ideal([f(x1+x2, x1-x2, 0), x3])
    G = curve_graph(I)
    @test number_of_connected_components(G) == 1

    ##########################
    ####### Edges cases ######

    R, (x1, x2, x3) = polynomial_ring(QQ, [:x1, :x2, :x3])

    # Zero Ideal
    I_zero = AlgebraicSolving.Ideal([R(0)])
    @test_throws Exception curve_graph(I_zero)

    # Empty Real Locus (x^2 + y^2 + 1 = 0)
    I_empty = AlgebraicSolving.Ideal([x1^2 + x2^2 + 1, x3])
    G_empty = curve_graph(I_empty)
    @test number_of_connected_components(G_empty) == 0
    @test length(G_empty.vertices) == 0

    # Simple (unbounded) straight line
    I_line = AlgebraicSolving.Ideal([x1 - x2 + 1, x3])
    G_line = curve_graph(I_line)
    @test number_of_connected_components(G_line) == 1

    # Output format flag (outf=false) testing exact arithmetic coordinates
    G_exact = curve_graph(I_line, outf=false)
    if length(G_exact.vertices) > 0
        # Check that the coordinates are exact fraction field elements, not Float64
        @test !(G_exact.vertices[1][1] isa Float64)
    end

    ##########################
    ####### Arrangements #####

    R, (x1, x2, x3) = polynomial_ring(QQ, [:x1, :x2, :x3])

    # Invalid Input (no curve)
    @test_throws AssertionError curve_arrangement_graph(typeof(AlgebraicSolving.Ideal([x1]))[])

    # Trivial Arrangement (Single Curve)
    I_line = AlgebraicSolving.Ideal([x2 - x3, x1])
    G_arr_single = curve_arrangement_graph([I_line])
    @test number_of_connected_components(G_arr_single) == 1

    # Two Disjoint Curves (Parallel Lines)
    I_line1 = AlgebraicSolving.Ideal([x2 - x3, x1])
    I_line2 = AlgebraicSolving.Ideal([x2 - x3 - 1, x1])
    G_arr_disjoint = curve_arrangement_graph([I_line1, I_line2])
    # They don't intersect, so they should remain 2 separate components
    @test number_of_connected_components(G_arr_disjoint) == 2

    # Two Intersecting Circles
    I_circle1 = AlgebraicSolving.Ideal([x2^2 + x3^2 - 1, x1])
    I_circle2 = AlgebraicSolving.Ideal([(x2-1)^2 + (x3-1)^2 - 2, x1])
    G_arr_intersect = curve_arrangement_graph([I_circle1, I_circle2])
    @test number_of_connected_components(G_arr_intersect) == 1

    # SKEW CURVES: Apparent Singularities
    # Their 2D projections cross at (0,0), but they DO NOT intersect in 3D.
    I_skew1 = AlgebraicSolving.Ideal([x2 - x3, x1 - 1])
    I_skew2 = AlgebraicSolving.Ideal([x2 + x3, x1 + 1])
    G_arr_skew = curve_arrangement_graph([I_skew1, I_skew2])
    @test number_of_connected_components(G_arr_skew) == 2
end