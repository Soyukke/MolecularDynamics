export LennardJones

"""
Potential model
"""
abstract type AbstractPotential end

function potential(p::T) where T <: AbstractPotential
    return p.f
end

struct LennardJones <: AbstractPotential
    f::Function
end

function LennardJones(;ε=one(Float64), σ=one(Float64))
    lj(x::T, y::T, z::T) where T = lennardjonespotential(x, y, z; ε, σ)
    return LennardJones(lj)
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
