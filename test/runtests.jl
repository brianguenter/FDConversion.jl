
import Symbolics
import FastDifferentiation as FD
using FDConversion
using Random
using Test


include("sphericalharmonics.jl")

@testset "conversion from Symbolics to FD" begin
    Symbolics.@variables x y

    symbolics_expr = x^2 + y * (x^2)
    dag, tmp = to_fd(symbolics_expr)
    vars = collect(values(tmp))
    fdx, fdy = FD.value(vars[1]) == :x ? (vars[1], vars[2]) : (vars[2], vars[1]) #need to find the variables since they can be in any order

    correct_fun = FD.make_function([dag], [fdx, fdy])


    #verify that all the node expressions exist in the dag. Can't rely on them being in a particular order because Symbolics can
    #arbitrarily choose how to reorder trees.
    num_tests = 100
    rng = Random.Xoshiro(8392)
    for _ in 1:num_tests
        (xval, yval) = rand(rng, 2)
        FDval = correct_fun([xval, yval])[1]
        Syval = Symbolics.substitute(symbolics_expr, Dict([(x, xval), (y, yval)]))

        @test isapprox(FDval, Syval.val)

    end
end


@testset "conversion from FD to Symbolics" begin
    order = 8
    FD.@variables x y z

    FD_funcs = FD.Node.(SHFunctions(order, x, y, z))
    Sym_funcs, variables = to_symbolics(FD_funcs)

    sx, sy, sz = map(p -> variables[p], [x, y, z])

    FD_eval = FD.make_function(FD_funcs, [x, y, z])
    rng = Random.Xoshiro(8392)
    num_tests = 100
    for _ in 1:num_tests
        tx, ty, tz = rand(rng, BigFloat, 3)
        subs = Dict([sx => tx, sy => ty, sz => tz])
        res = Symbolics.substitute.(Sym_funcs, Ref(subs))

        FD_res = FD_eval([tx, ty, tz])

        @test isapprox(FD_res, res, atol=1e-12)
    end
end


