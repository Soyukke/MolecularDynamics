export forcefunc

"""
Force model
"""
abstract type AbstractForce end

"""
Calculate force function by potential function
"""
function forcefunc(p::T) where T <: AbstractPotential
    pot = potential(p)
    return force(x, y, z) = gradient(x, y, z) do x₁, y₁, z₁
        pot(x₁, y₁, z₁)
    end
end
