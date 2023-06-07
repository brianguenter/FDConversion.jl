module FDConversion

using StaticArrays
import Symbolics
import SymbolicUtils
import FastDifferentiation as FD
using FastDifferentiation.AutomaticDifferentiation

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

function to_FD(x::Real, cache::IdDict=IdDict(), substitutions::Union{IdDict,Nothing}=nothing)
    return _to_FD(x, cache, substitutions)
end
export to_FD


#WARNING!!!!!!! TODO. *,+ simplification code relies on Node(0) have node_value 0 as Int64. Need to make sure that expr_to_dag unwraps numbers or the simplification code won't work.
function _to_FD(symx, cache::IdDict, substitutions::Union{IdDict,Nothing}) #cache is an IdDict, to make clear that hashing into the cache Dict is  based on objectid, i.e., using === rather than ==.
    # Substitutions are done on a Node graph, not a SymbolicsUtils.Sym graph. As a consequence the values
    # in substitutions Dict are Node not Sym type. cache has keys (op,args...) where op is generally a function type but sometimes a Sym, 
    # and args are all Node types.

    #need to extract the SymbolicUtils tree from symx

    if isa(symx, Num) #substitutions are stored as SymbolicUtils.Symx so extract the underlying Symx value
        symx = symx.val
    end

    if substitutions !== nothing
        tmpsub = get(substitutions, symx, nothing)
        if tmpsub !== nothing
            return substitutions[symx] #substitute Node object for symbolic object created in differentiation
        end
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


        #Taking Ref(nothing) causes the broadcasting to screw up for some reason. need two cases on for subsitutions === nothing and one for it being an IdDict.
        if substitutions === nothing
            args = _to_FD.(symargs, Ref(cache), substitutions)
        else
            args = _to_FD.(symargs, Ref(cache), Ref(substitutions))
        end

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

function P(l, m, z)
    if l == 0 && m == 0
        return 1.0
    elseif l == m
        return (1 - 2m) * P(m - 1, m - 1, z)
    elseif l == m + 1
        return (2m + 1) * z * P(m, m, z)
    else
        return ((2l - 1) / (l - m) * z * P(l - 1, m, z) - (l + m - 1) / (l - m) * P(l - 2, m, z))
    end
end


function S(m, x, y)
    if m == 0
        return 0
    else
        return x * C(m - 1, x, y) - y * S(m - 1, x, y)
    end
end


function C(m, x, y)
    if m == 0
        return 1
    else
        return x * S(m - 1, x, y) + y * C(m - 1, x, y)
    end
end

function factorial_approximation(x)
    local n1 = x
    sqrt(2 * π * n1) * (n1 / ℯ * sqrt(n1 * sinh(1 / n1) + 1 / (810 * n1^6)))^n1
end


function compare_factorial_approximation()
    for n in 1:30
        println("n $n relative error $((factorial(big(n))-factorial_approximation(n))/factorial(big(n)))")
    end
end


function N(l, m)
    @assert m >= 0
    if m == 0
        return sqrt((2l + 1 / (4π)))
    else
        # return sqrt((2l+1)/2π * factorial(big(l-m))/factorial(big(l+m)))
        #use factorial_approximation instead of factorial because the latter does not use Stirlings approximation for large n. Get error for n > 2 unless using BigInt but if use BigInt get lots of rational numbers in symbolic result.
        return sqrt((2l + 1) / 2π * factorial_approximation(l - m) / factorial_approximation(l + m))
    end
end


"""l is the order of the spherical harmonic"""
function Y(l, m, x, y, z)
    @assert l >= 0
    @assert abs(m) <= l
    if m < 0
        return N(l, abs(m)) * P(l, abs(m), z) * S(abs(m), x, y)
    else
        return N(l, m) * P(l, m, z) * C(m, x, y)
    end
end


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

    from_dag = to_symbolics.(FD_funcs)
    for i in 1:1_000
        tx, ty, tz = rand(BigFloat, 3)
        subs = Dict([sx => tx, sy => ty, sz => tz])
        @assert isapprox(FD_eval([tx, ty, tz]), Symbolics.substitute.(Sym_funcs, Ref(subs)), atol=1e-12)
    end
end
export test

end # module FSDConvert
