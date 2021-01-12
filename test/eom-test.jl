using MolecularDynamics
using LinearAlgebra
using Test

function testfractionalcoords01()
    cell = Float64[10, 10, 10]
    A = Matrix{Float64}(Diagonal(cell))
    x = Float64[3, 4, 5]
    x̂ = MolecularDynamics.fractionalcoords(A, x)
    return round.(x̂, digits=3) == [0.3, 0.4, 0.5]
end

function testfractionalcoords02()
    cell = Float64[10, 10, 10]
    A = Matrix{Float64}(Diagonal(cell))
    x = Float64[3 1; 4 2; 5 3;]
    x̂ = MolecularDynamics.fractionalcoords(A, x)
    return round.(x̂, digits=3) == [0.3 0.1; 0.4 0.2; 0.5 0.3]
end

function testcartesiancoords01()
    cell = Float64[10, 10, 10]
    A = Matrix{Float64}(Diagonal(cell))
    x̂ = Float64[0.3, 0.4, 0.5]
    x = MolecularDynamics.cartesiancoords(A, x̂)
    return round.(x, digits=3) == [3, 4, 5]
end

function testcartesiancoords02()
    cell = Float64[10, 10, 10]
    A = Matrix{Float64}(Diagonal(cell))
    x̂ = Float64[0.3 0.1; 0.4 0.2; 0.5 0.3]
    x = MolecularDynamics.cartesiancoords(A, x̂)
    return round.(x, digits=3) == [3 1; 4 2; 5 3]
end

@testset "eom.jl" begin
    @test testfractionalcoords01()
    @test testfractionalcoords02()
    @test testcartesiancoords01()
    @test testcartesiancoords02()
end