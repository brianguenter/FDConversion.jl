
function test1()
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
        println(FD_res - res)
        @assert isapprox(FD_res, res, atol=1e-12)
    end
end
export test1