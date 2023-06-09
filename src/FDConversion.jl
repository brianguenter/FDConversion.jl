module FDConversion

using StaticArrays
import Symbolics
import SymbolicUtils
import FastDifferentiation as FD
using FastDifferentiation.AutomaticDifferentiation
import Random

"""converts from Node to Symbolics expression"""
function to_symbolics(a::T, cache::IdDict=IdDict()) where {T<:FD.Node}
    tmp = get(cache, a, nothing)
    if tmp !== nothing
        return tmp
    else
        if FD.arity(a) === 0
            if FD.is_constant(a)
                cache[a] = Num(FD.value(a))
            elseif FD.is_variable(a)
                cache[a] = Symbolics.variable(FD.value(a))
            else
                throw(ErrorException("Node with 0 children was neither constant nor variable. This should never happen."))
            end
        else
            if FD.arity(a) === 1
                cache[a] = a.node_value(to_symbolics(a.children[1], cache))
            else
                cache[a] = foldl(a.node_value, to_symbolics.(a.children, Ref(cache)))
            end
        end

        return cache[a]
    end
end
export to_symbolics

to_symbolics(a::AbstractVector{T}, cache::IdDict=IdDict()) where {T<:FD.Node} = to_symbolics.(a, Ref(cache))

to_FD(x::FD.AutomaticDifferentiation.NoDeriv, cache, substitions) = FD.Node(NaN) #when taking the derivative with respect to the first element of 1.0*x Symbolics.derivative will return Symbolics.NoDeriv. These derivative values will never be used (or should never be used) in my derivative computation so set to NaN so error will show up if this ever happens.

function to_FD(x::Real)
    cache = IdDict()
    substitutions = IdDict()
    return _to_FD(x, cache, substitutions)
end
export to_FD


function _to_FD(sym_node::Real, cache::IdDict, visited::IdDict)
    # Substitutions are done on a Node graph, not a SymbolicsUtils.Sym graph. As a consequence the values
    # in substitutions Dict are Node not Sym type. cache has keys (op,args...) where op is generally a function type but sometimes a Sym, 
    # and args are all Node types.

    symx = sym_node.val

    tmpsub = get(visited, symx, nothing)
    if tmpsub !== nothing
        return visited[symx] #substitute Node object for symbolic object
    end

    tmp = get(cache, symx, nothing)

    if tmp !== nothing
        return tmp
    elseif !SymbolicUtils.istree(symx)

        tmpnode = FD.Node(symx)
        cache[symx] = tmpnode

        return tmpnode
    else
        numargs = length(SymbolicUtils.arguments(symx))
        symargs = MVector{numargs}(SymbolicUtils.arguments(symx))

        args = _to_FD.(symargs, Ref(cache), Ref(visited))

        key = (SymbolicUtils.operation(symx), args...)
        tmp = get(cache, key, nothing)

        if tmp !== nothing
            return tmp
        else
            tmpnode = FD.Node(SymbolicUtils.operation(symx), args)
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

        @assert isapprox(FD_res, res, atol=1e-12) "error $(FD_res - res)"
    end
end
export test

end # module FSDConvert
