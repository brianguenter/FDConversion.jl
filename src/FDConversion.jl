module FDConversion

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

end # module FSDConvert
