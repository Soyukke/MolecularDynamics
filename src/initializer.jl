abstract type AbstractInitializer end
abstract type DISTRIBUTE_MAXWELL <: AbstractInitializer end
abstract type DISTRIBUTE_UNIFORM <: AbstractInitializer end

"""
init_vector

Initialize vector by distributions.
"""
function init_vector(shape::Tuple, inittype::Type{T}) where T <: AbstractInitializer
    #= init vector value =#
    if inittype == DISTRIBUTE_MAXWELL
        return randn(shape...)
    elseif inittype == DISTRIBUTE_UNIFORM
        return rand(shape...)
    end
end