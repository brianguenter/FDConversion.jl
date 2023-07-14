
module ConversionTests

using TestItems

@testitem "to_FD" begin

    import Symbolics
    import FastDifferentiation as FD
    using Random

    Symbolics.@variables x y

    symbolics_expr = x^2 + y * (x^2)
    dag = to_fd(symbolics_expr)
    vars = FD.variables(dag)
    fdx, fdy = FD.value(vars[1]) == :x ? (vars[1], vars[2]) : (vars[2], vars[1]) #need to find the variables since they can be in any order

    correct_fun = FD.make_function([dag], [fdx, fdy])


    #verify that all the node expressions exist in the dag. Can't rely on them being in a particular order because Symbolics can
    #arbitrarily choose how to reorder trees.
    num_tests = 1_000
    rng = Random.Xoshiro(8392)
    for _ in 1:num_tests
        (xval, yval) = rand(rng, 2)
        FDval = correct_fun([xval, yval])[1]
        Syval = Symbolics.substitute(symbolics_expr, Dict([(x, xval), (y, yval)]))

        @test isapprox(FDval, Syval.val)

    end
end


@testitem "to_symbolics" begin #test conversion from FD to Symbolics
    import FastDifferentiation as FD
    import Symbolics
    import Random

    order = 8
    FD.@variables x y z

    FD_funcs = FD.Node.(SHFunctions(order, x, y, z))
    Sym_funcs, variables = to_symbolics(FD_funcs)

    sx, sy, sz = map(p -> variables[p], [x, y, z])

    FD_eval = FD.make_function(FD_funcs, [x, y, z])
    rng = Random.Xoshiro(8392)
    for _ in 1:1_0
        tx, ty, tz = rand(rng, BigFloat, 3)
        subs = Dict([sx => tx, sy => ty, sz => tz])
        res = Symbolics.substitute.(Sym_funcs, Ref(subs))

        FD_res = FD_eval([tx, ty, tz])

        @assert isapprox(FD_res, res, atol=1e-12)
    end
end