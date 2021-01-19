export LennardJones

"""
Potential model
"""
abstract type AbstractPotential end

function potential(p::T) where T <: AbstractPotential
    return p.f
end

"""
General potential function struct
"""
struct PotentialFunction <: AbstractPotential
    f::Function
end

struct LennardJones <: AbstractPotential
    f::Function
end

"""
LennardJones(;ε=one(Float64), σ=one(Float64))

`ε` and `σ` are parameters
"""
function LennardJones(;ε=one(Float64), σ=one(Float64))
    lj(x::T, y::T, z::T) where T = lennardjonespotential(x, y, z; ε, σ)
    return PotentialFunction(lj)
end

"""
potential
user setting potential function
not use potential for calculation of step
"""
function lennardjonespotential(x::T, y::T, z::T; ε=one(T), σ=one(T)) where T
    r = sqrt(x^2 + y^2 + z^2)
    x = σ/r
    b = x^6
    a = b^2
    return 4ε*(a - b)
end


"""
Atomic pair force function
For create potential settings

TODO For use sum potential functions and calculate force function
"""
struct PotentialMap
    # {atomic type, atomic type} => ForcuFunction
    dict::Dict{Pair{Integer, Integer}, PotentialFunction}
end

"""
get potential function by atomic type pair
i and j == j and i
"""
function (m::PotentialMap)(i::Integer, j::Integer)
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