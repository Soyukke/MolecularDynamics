export forcefunc, generalforcefunction, ForceMap

"""
Force model
"""
abstract type AbstractForce end

struct ForceFunction <: AbstractForce
    f::Function
end

function generalforcefunction(f::Function)
    return ForceFunction(f)
end

"""
Calculate force function by potential function
"""
function forcefunc(p::T) where T <: AbstractPotential
    pot = potential(p)
    return force(x, y, z) = gradient(x, y, z) do x₁, y₁, z₁
        - pot(x₁, y₁, z₁)
    end
end

"""
Atomic pair force function
"""
struct ForceMap
    # {atomic type, atomic type} => ForcuFunction
    dict::Dict{Pair{Integer, Integer}, ForceFunction}
end

"""
get force function by atomic type pair
i and j == j and i
"""
function (m::ForceMap)(i::Integer, j::Integer)
    key = 
    if haskey(m.dict, Pair(i, j))
        Pair(i, j)
    elseif haskey(m.dict, Pair(j, i))
        Pair(j, i)
    else
        nothing
    end
    return !isnothing(key) ? m.dict[key].f : nothing
end