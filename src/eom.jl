using Printf
using LinearAlgebra
using Base.Iterators
using Base.Threads
using Combinatorics
using Random
import Base.copy
const INIT_MAXWELL = "maxwell"
const INIT_UNIFORM = "random"

function printarray(x)
    show(stdout, "text/plain", x)
    @printf "\n"
end

mutable struct Geometry
    nstep::Int64
    Δt::Float64 # time step
    cell::Array{Float64, 1} # x, y, z
    natom::Int64 # number of particles
    vnatom::Array{Int64, 1} # atom of type
    ntype::Int64
    vtype::Array{Int64, 1} # type
    vm::Array{Float64, 2} # mass
    vr::Array{Float64, 2} # position
    vv::Array{Float64, 2} # velocity
    # force
    R1::Float64
    R2::Float64
    R::Float64
    Life::Float64
    vflist::Array{Array{Int64, 1}, 1} # 粒子のインデックスと距離を保存する
    vf0::Array{Float64, 2}
    vf::Array{Float64, 2}
end

"""
copy(g::Geometry)

copy `g` Geometry
"""
function copy(g::Geometry)
    return Geometry(
        g.nstep,
        g.Δt,
        g.cell,
        g.natom,
        g.vnatom,
        g.ntype,
        g.vtype,
        copy(g.vm),
        copy(g.vr),
        copy(g.vv),
        g.R1,
        g.R2,
        g.R,
        g.Life,
        copy(g.vflist),
        copy(g.vf0),
        copy(g.vf),
    )
end

"""
potential
user setting potential function
not use potential for calculation of step
"""
function potential(r)
    ε = 1
    σ = 1
    x = σ/r
    b = x^6
    a = b^2
    return 4ε*(a - b)
end

"""
force
force function
calculate force by potential
"""
function force(r)
    # norm of force
    # vector is \vec{e_r} = x/sqrt(x^2 + y^2 + z^2)e_x + y/....
    # vF = ∇V
    ε = 1
    σ = 1
    x = σ/r
    # diff x = -σ/r^2 = - x^2 / σ
    b = -6*x^7 / σ
    a = -12x^13 / σ
    return -4ε*(a - b)
end

"""
init_vector

Initialize vector by distributions.
"""
function init_vector(shape::Tuple, init_type=INIT_MAXWELL)
    #= init vector value =#
    if init_type == INIT_MAXWELL
        return randn(shape...)
    elseif init_type == INIT_UNIFORM
        return rand(shape...)
    end
end

function sample2hist(vx::Array{Float64, 1}, dt::Float64)
    lvx = sort(vx)
    max_vx = max(vx...)
    min_vx = min(vx...)
    vt = min_vx:dt:max_vx
    vf = zeros(length(vt))
    vf = map(t -> length(lvx[t .< lvx .< t+dt]), vt)
    plot(vt, abs.(vt).*vf)
    savefig("vmaxwell.png")
end

"""
make_force_list
geometryのforcelistを再計算する
一定距離以内の原子をmemoryすることで計算量を減らす
"""
function make_force_list(geometry::Geometry)
    natom = geometry.natom
    vr = geometry.vr
    # いったん初期化
    geometry.vflist = Array{Array{Int64, 1}, 1}[]
    vflist = geometry.vflist
    for i in 1:natom push!(vflist, []) end
    R = geometry.R
    # 原子ペアの全組み合わせ
    pairs = combinations(1:natom, 2)
    # TODO parallelable
    for (i, j) in pairs
        r = norm(vr[:, i] - vr[:, j])
        if r < R
            push!(vflist[i], j)
        end
    end
end

"""
make_geometry

make example geometry
"""
function make_geometry()
    # ジオメトリ作成
    nstep = 1000
    Δt = 1e-2
    cell = [10, 10, 10]
    atomnumbers = [1, 6]
    # 各原子数
    vnatom = [6, 12]
    natom = sum(vnatom)
    ntype = 2
    # [1, 1, 1, 1,  1,  1,  1,  1,  ... 2,  2,  2,  2,  2]
    vtype = vcat(map(x -> x[1]*ones(Int64, x[2]), enumerate(vnatom))...)
    vm = ones(natom)'
    vr = 10*init_vector((3, natom), INIT_UNIFORM)
    vv = zeros(3, natom)
    # force
    R1 = 4
    R2 = 4
    R = 8
    L = 0
    vlist = Array{Array{Int64, 1}, 1}[]
    for i in 1:natom push!(vlist, []) end
    vf0 = zeros(3, natom)
    vf = zeros(3, natom)
    return Geometry(nstep, Δt, cell, natom, atomnumbers, ntype, vtype, vm, vr, vv, R1, R2, R, L, vlist, vf0, vf)
end

"""
step

calculate 1 step
"""
function step(geometry::Geometry)
    natom = geometry.natom
    Δt = geometry.Δt
    vm = geometry.vm
    vr = geometry.vr
    vv = geometry.vv
    vflist = geometry.vflist
    vf0 = geometry.vf0
    vf = geometry.vf
    vr .= vr .+ vv*Δt + (Δt^2 ./ 2vm .* vf)
    vf0 .= copy(vf)
    vf .= zeros(3, natom)

    for (i, vj) in enumerate(vflist)
        for j in vj
            r_12 = vr[:, i] - vr[:, j]
            r = norm(r_12)
            nr_12 = r_12 / r
            vf[:, i] += force(r)*nr_12
            vf[:, j] -= force(r)*nr_12
        end
    end

    # calculate vv
    vv .= vv + (Δt^2 ./ 2vm .* (vf0 + vf))
end

function main()
    geom = make_geometry()
    geoms = Geometry[geom]
    make_force_list(copy(geom))
    for i in 1:1000
        step(geom)
        push!(geoms, copy(geom))
        geom.Life += max(mapslices(norm, geom.vr, dims=1)...) * geom.Δt
        # @info geom.Life, geom.R2
        if geom.Life > geom.R2
            make_force_list(geom)
            geom.Life = 0
        end
    end
    return geoms
end
