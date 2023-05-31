
@testitem "all_nodes" begin
    using FastSymbolicDifferentiation
    import Symbolics

    Symbolics.@variables x y
    expr_to_dag(x^2 + y * (x^2), cache), x, y
    cache = IdDict()
    dag, x, y = simple_dag(cache)

    correct = expr_to_dag.([((x^2) + (y * (x^2))), (x^2), x, 2, (y * (x^2)), y], Ref(cache))
    tmp = all_nodes(dag)

    #verify that all the node expressions exist in the dag. Can't rely on them being in a particular order because Symbolics can
    #arbitrarily choose how to reorder trees.
    for expr in correct
        @test in(expr, tmp)
    end
end


@testitem "conversion from graph of FastSymbolicDifferentiation.Node to Symbolics expression" begin
    using FastSymbolicDifferentiation
    import Symbolics

    order = 7
    @variables x y z

    derivs = Symbolics.jacobian(SHFunctions(order, x, y, z), [x, y, z]; simplify=true)
    # show(@time SHDerivatives(order,x,y,z))
    tmp = expr_to_dag.(derivs)
    # show(@time expr_to_dag.(derivs))
    from_dag = to_symbolics.(tmp)
    subs = Dict([x => rand(), y => rand(), z => rand()])
    @test isapprox(map(xx -> xx.val, Symbolics.substitute.(derivs, Ref(subs))), map(xx -> xx.val, Symbolics.substitute.(from_dag, Ref(subs))), atol=1e-12)
end