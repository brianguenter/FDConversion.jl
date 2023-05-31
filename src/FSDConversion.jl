module FSDConvert

using FastSymbolicDifferentiation
using Symbolics

"""converts from dag to Symbolics expression"""
function to_symbolics(a::Node)
    if arity(a) === 0
        return Num(value(a)) #convert everything to Num type. This will wrap types like Int64,Float64, etc., but will not double wrap nodes that are Num types already.
    else
        if arity(a) === 1
            return a.node_value(to_symbolics(a.children[1]))
        else
            return foldl(a.node_value, to_symbolics.(a.children))
        end
    end
end
export to_symbolics

expr_to_dag(x::AutomaticDifferentiation.NoDeriv, cache, substitions) = Node(NaN) #when taking the derivative with respect to the first element of 1.0*x Symbolics.derivative will return Symbolics.NoDeriv. These derivative values will never be used (or should never be used) in my derivative computation so set to NaN so error will show up if this ever happens.

function expr_to_dag(x::Real, cache::IdDict=IdDict(), substitutions::Union{IdDict,Nothing}=nothing)
    return _expr_to_dag(x, cache, substitutions)
end
export expr_to_dag


#WARNING!!!!!!! TODO. *,+ simplification code relies on Node(0) have node_value 0 as Int64. Need to make sure that expr_to_dag unwraps numbers or the simplification code won't work.
function _expr_to_dag(symx, cache::IdDict, substitutions::Union{IdDict,Nothing}) #cache is an IdDict, to make clear that hashing into the cache Dict is  based on objectid, i.e., using === rather than ==.
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

        tmpnode = Node(symx)
        cache[symx] = tmpnode

        return tmpnode
    else
        numargs = length(SymbolicUtils.arguments(symx))
        symargs = MVector{numargs}(SymbolicUtils.arguments(symx))


        #Taking Ref(nothing) causes the broadcasting to screw up for some reason. need two cases on for subsitutions === nothing and one for it being an IdDict.
        if substitutions === nothing
            args = _expr_to_dag.(symargs, Ref(cache), substitutions)
        else
            args = _expr_to_dag.(symargs, Ref(cache), Ref(substitutions))
        end

        key = (SymbolicUtils.operation(symx), args...)
        tmp = get(cache, key, nothing)

        if tmp !== nothing
            return tmp
        else
            tmpnode = Node(SymbolicUtils.operation(symx), args)
            cache[key] = tmpnode

            return tmpnode
        end
    end
end

end # module FSDConvert
