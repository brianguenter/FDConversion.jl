
module ConversionTests

using TestItems

@testitem "all_nodes" begin
    # Testing for equivalence of expression graphs is inherently difficult. 
    # Symbolics and FastDifferentiation apply different ordering rules and simplifications at expression construction time
    # and these rules may change over time. Example:

    # Symbolics reorders terms in * expressions

    # Symbolics.@variables x y
    # julia> x*y
    # x*y
    # julia> y*x
    # x*y

    # FastDifferentiation does not reorder terms as of v0.2.0

    # julia> FastDifferentiation.@variables x1 y1
    # julia> x1*y1
    # (x1 * y1)
    # julia> y1*x1
    # (y1 * x1)

    # Most reliable way to test for equivalence is to evaluate the function at many points and test for approximate equality of result.

    import FastDifferentiation as FD
    import Symbolics

    Symbolics.@variables x y

    symbolics_expr = x^2 + y * (x^2)
    dag = to_FD(symbolics_expr)
    fdx, fdy = FD.variables(dag)

    correct = fdx^2 + fdy * (fdx^2)
    correct_fun = FD.make_function([correct], [fdx, fdy])

    #verify that all the node expressions exist in the dag. Can't rely on them being in a particular order because Symbolics can
    #arbitrarily choose how to reorder trees.
    num_tests = 1_000
    rng = Random.Xoshiro(8392)
    for _ in num_tests
        for (xval, yval) in rand(rng, 2)
            FDval = correct_fun(xval, yval)
            Syval = Symbolics.substitute(symbolics_expr, Dict([(x, xval), (y, yval)]))
            @test isapprox(FDval, Syval)
        end
    end
end


@testitem "to_symbolics" begin #test conversion from FD to Symbolics
    import FastDifferentiation as FD
    import Symbolics
    import Random

    order = 8
    FD.@variables x y z
    Symbolics.@variables sx, sy, sz

    FD_funcs = FD.Node.(SHFunctions(order, x, y, z))
    Sym_funcs = SHFunctions(order, sx, sy, sz)

    FD_eval = FD.make_function(FD_funcs, [x, y, z])
    rng = Random.Xoshiro(8392)
    for _ in 1:1_0
        tx, ty, tz = rand(rng, BigFloat, 3)
        subs = Dict([sx => tx, sy => ty, sz => tz])
        res = Symbolics.substitute.(Sym_funcs, Ref(subs))

        FD_res = FD_eval([tx, ty, tz])

        @test isapprox(FD_res, res, atol=1e-12)
    end
end

end #module