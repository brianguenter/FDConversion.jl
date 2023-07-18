module FDConversion

import Symbolics
import SymbolicUtils
import FastDifferentiation as FD
using FastDifferentiation.AutomaticDifferentiation
import Random

#=
try implementing these
    https://github.com/JuliaSymbolics/SymbolicUtils.jl/blob/master/src/interface.jl
=#

"""converts from Node to Symbolics expression"""
function _to_symbolics!(a::T, cache::IdDict, variable_map::IdDict) where {T<:FD.Node}
    tmp = get(cache, a, nothing)
    if tmp !== nothing
        return tmp
    else
        if FD.arity(a) === 0
            if FD.is_constant(a)
                cache[a] = Symbolics.Num(FD.value(a))
            elseif FD.is_variable(a)
                tmp = Symbolics.variable(FD.value(a))
                cache[a] = tmp
                variable_map[a] = tmp
            else
                throw(ErrorException("Node with 0 children was neither constant nor variable. This should never happen."))
            end
        else
            if FD.arity(a) === 1
                cache[a] = a.node_value(_to_symbolics!(a.children[1], cache, variable_map))
            else
                cache[a] = foldl(a.node_value, _to_symbolics!.(a.children, Ref(cache), Ref(variable_map)))
            end
        end
    end
end

function to_symbolics(a::T) where {T<:FD.Node}
    cache = IdDict()
    variable_map = IdDict()
    _to_symbolics!(a, cache, variable_map)
    return cache[a], variable_map
end
export to_symbolics

function to_symbolics(a::AbstractArray{T}) where {T<:FD.Node}
    cache::IdDict = IdDict()
    variable_map::IdDict = IdDict()

    _to_symbolics!.(a, Ref(cache), Ref(variable_map))
    return map(x -> cache[x], a), variable_map
end


function to_fd(x::Real, cache=IdDict(), substitutions=IdDict())
    result = _to_FD!(x, cache, substitutions)
    syms = collect(filter(x -> SymbolicUtils.issym(x), keys(cache)))
    sym_map = Dict(zip(syms, map(x -> cache[x], syms)))

    return result, sym_map
end
export to_fd

function to_fd(x::AbstractArray{<:Real})
    cache = IdDict()
    substitutions = IdDict()
    result = _to_FD!.(x, Ref(cache), Ref(substitutions))
    syms = collect(filter(x -> SymbolicUtils.issym(x), keys(cache)))
    sym_map = Dict(zip(syms, map(x -> cache[x], syms))) #map from symbolics variables to FD variables
    return result, sym_map
end

function _to_FD!(sym_node, cache::IdDict, visited::IdDict)
    # Substitutions are done on a Node graph, not a SymbolicsUtils.Sym graph. As a consequence the values
    # in substitutions Dict are Node not Sym type. cache has keys (op,args...) where op is generally a function type but sometimes a Sym, 
    # and args are all Node types.

    @assert !(typeof(sym_node) <: Symbolics.Arr) "Differentiation of expressions involving arrays and array variables is not yet supported."

    symx = isa(sym_node, Symbolics.Num) ? sym_node.val : sym_node
    @assert typeof(symx) != Symbolics.Num


    tmpsub = get(visited, symx, nothing)
    if tmpsub !== nothing
        return visited[symx] #substitute Node object for symbolic object
    end

    tmp = get(cache, symx, nothing)

    if tmp !== nothing
        return tmp
    elseif !SymbolicUtils.istree(symx)
        if SymbolicUtils.issym(symx)
            tmpnode = FD.Node(Symbol(symx))
        else #must be a number of some kind
            tmpnode = FD.Node(symx)
        end

        cache[symx] = tmpnode

        return tmpnode
    else
        numargs = length(SymbolicUtils.arguments(symx))
        symargs = SymbolicUtils.arguments(symx)

        args = _to_FD!.(symargs, Ref(cache), Ref(visited))

        key = (SymbolicUtils.operation(symx), args...)
        tmp = get(cache, key, nothing)

        if tmp !== nothing
            return tmp
        else
            tmpnode = FD.Node(SymbolicUtils.operation(symx), args...)
            cache[key] = tmpnode

            return tmpnode
        end
    end
end


"""
Converts from `Symbolics` form to `FastDifferentiation` form and computes Jacobian with respect to `diff_variables`.
If `fast_differentiation=false` the result will be in Symbolics form. If `fast_differentiation=true` 
then the result will be a two tuple. The first tuple entry will be `function` converted to `FastDifferentiation` form. 
The second tuple term will be `diff_variables` converted to `FastDifferentiation` form.
These two values can then be used to make an efficient executable using `make_function`."""
function fd_jacobian(symbolics_function::AbstractArray{Symbolics.Num}, differentiation_variables::AbstractVector{Symbolics.Num}, fast_differentiation=false)
    fd_func = to_fd(func)
    tmp = jacobian(fd_func)
end


function fd_sparse_jacobian()
end
function fd_hessian()
end
function fd_sparse_hessian()
end
#etc. for Jv Jᵀv Hv

function test()
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

        @assert isapprox(FDval, Syval.val)

    end
end
export test
# export FastDifferentiation.make_function
end # module FSDConvert

