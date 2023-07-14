module FDConversion

using StaticArrays
import Symbolics
import SymbolicUtils
import FastDifferentiation as FD
using FastDifferentiation.AutomaticDifferentiation
import Random

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
export _to_symbolics!

function to_symbolics(a::T) where {T<:FD.Node}
    cache = IdDict()
    variable_map = IdDict()
    _to_symbolics!(a, cache, variable_map)
    return cache[a], variable_map
end

function to_symbolics(a::AbstractVector{T}) where {T<:FD.Node}
    cache::IdDict = IdDict()
    variable_map::IdDict = IdDict()

    _to_symbolics!.(a, Ref(cache), Ref(variable_map))
    return map(x -> cache[x], a), variable_map
end

to_fd(x::FD.AutomaticDifferentiation.NoDeriv, cache, substitions) = FD.Node(NaN) #when taking the derivative with respect to the first element of 1.0*x Symbolics.derivative will return Symbolics.NoDeriv. These derivative values will never be used (or should never be used) in my derivative computation so set to NaN so error will show up if this ever happens.


function to_fd(x::Real)
    cache = IdDict()
    substitutions = IdDict()
    return _to_FD(x, cache, substitutions)
end
export to_fd


function _to_FD(sym_node, cache::IdDict, visited::IdDict)
    # Substitutions are done on a Node graph, not a SymbolicsUtils.Sym graph. As a consequence the values
    # in substitutions Dict are Node not Sym type. cache has keys (op,args...) where op is generally a function type but sometimes a Sym, 
    # and args are all Node types.

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

        args = _to_FD.(symargs, Ref(cache), Ref(visited))

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

function P(::Type{T}, l, m, z::T) where {T}
    if l == 0 && m == 0
        return T(1)
    elseif l == m
        return (1 - 2m) * P(T, m - 1, m - 1, z)
    elseif l == m + 1
        return (2m + 1) * z * P(T, m, m, z)
    else
        return ((2l - 1) / (l - m) * z * P(T, l - 1, m, z) - (l + m - 1) / (l - m) * P(T, l - 2, m, z))
    end
end


function S(::Type{T}, m, x::T, y::T) where {T}
    if m == 0
        return T(0)
    else
        return x * C(T, m - 1, x, y) - y * S(T, m - 1, x, y)
    end
end


function C(::Type{T}, m, x::T, y::T) where {T}
    if m == 0
        return T(1)
    else
        return x * S(T, m - 1, x, y) + y * C(T, m - 1, x, y)
    end
end

function factorial_approximation(::Type{T}, x) where {T}
    local n1 = x
    sqrt(2 * T(π) * n1) * (n1 / T(ℯ) * sqrt(n1 * sinh(1 / T(n1)) + 1 / (810 * T(n1)^6)))^n1
end


function compare_factorial_approximation()
    for n in 1:30
        println("n $n relative error $((factorial(big(n))-factorial_approximation(BigFloat,n))/factorial(big(n)))")
    end
end


function N(::Type{T}, l, m) where {T}
    @assert m >= 0
    if m == 0
        return sqrt((2l + 1 / (4 * T(π))))
    else
        # return sqrt((2l+1)/2π * factorial(big(l-m))/factorial(big(l+m)))
        #use factorial_approximation instead of factorial because the latter does not use Stirlings approximation for large n. Get error for n > 2 unless using BigInt but if use BigInt get lots of rational numbers in symbolic result.
        return sqrt((2l + 1) / 2 * T(π) * factorial_approximation(T, l - m) / factorial_approximation(T, l + m))
    end
end


"""l is the order of the spherical harmonic"""
function Y(::Type{T}, l, m, x::T, y::T, z::T) where {T}
    @assert l >= 0
    @assert abs(m) <= l
    if m < 0
        return N(T, l, abs(m)) * P(T, l, abs(m), z) * S(T, abs(m), x, y)
    else
        return N(T, l, m) * P(T, l, m, z) * C(T, m, x, y)
    end
end

Y(l, m, x::T, y::T, z::T) where {T<:FD.Node} = Y(FD.Node, l, m, x, y, z)
Y(l, m, x::T, y::T, z::T) where {T<:Number} = Y(T, l, m, x, y, z)


function SHFunctions(max_l, x, y, z)
    @assert typeof(x) == typeof(y) == typeof(z)
    result = Vector(undef, 0)

    for l in 0:max_l-1
        for m in -l:l
            push!(result, Y(l, m, x, y, z))
        end
    end

    return result
end
export SHFunctions

function test()
    order = 8
    Symbolics.@variables sx, sy, sz

    Sym_funcs = SHFunctions(order, sx, sy, sz)
    FD_funcs = to_fd.(SymFuncs)
    fd_vars = FD.variables(FD_funcs)


    FD_eval = FD.make_function(FD_funcs, [x, y, z])
    rng = Random.Xoshiro(8392)
    for _ in 1:1_0
        tx, ty, tz = rand(rng, BigFloat, 3)
        subs = Dict([sx => tx, sy => ty, sz => tz])
        res = Symbolics.substitute.(Sym_funcs, Ref(subs))

        FD_res = FD_eval([tx, ty, tz])

        @assert isapprox(FD_res, res, atol=1e-12) "error $(FD_res - res)"
    end
end
export test

function test2()

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

        @assert isapprox(FDval, Syval.val)

    end
end
export test2

include("scratchpad.jl")
end # module FSDConvert
