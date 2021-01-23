using Printf
using LinearAlgebra
using Base.Iterators
using Base.Threads
using Combinatorics
using Random
import Base.copy

function printarray(x)
    show(stdout, "text/plain", x)
    @printf "\n"
end

mutable struct Geometry
    nstep::Int64
    step::Int64
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
    vflistexpand::Vector{Vector{Tuple{Int64, Tuple{Int64, Int64, Int64}}}} # 距離R以内の粒子のインデックスとイメージインデックスを保存する
    S::Matrix{Float64} # vflistexpandを更新してからの原子シフト
    vf0::Array{Float64, 2}
    vf::Array{Float64, 2}
    isperiodic::Vector{Bool}
    forcemap::ForceMap
end

"""
copy(g::Geometry)

copy `g` Geometry
"""
function copy(g::Geometry)
    return Geometry(
        g.nstep,
        g.step,
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
        copy(g.vflistexpand),
        copy(g.S),
        copy(g.vf0),
        copy(g.vf),
        copy(g.isperiodic),
        g.forcemap
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
make_force_list_expand


周期境界ありgeometryのforcelistを再計算する
一定距離以内の原子をmemoryすることで計算量を減らす

粒子番号だけでは情報が足りない
粒子番号, (IMAGE_X, IMAGE_Y, IMAGE_Z)
i, (+1, -1, +1)
Tuple{Int64, Int64, Int64} は 原子座標が周期境界で変更した場合に修正する
"""
function forcelistexpand(geometry::Geometry)
    # periodically
    isloopx, isloopy, isloopz = geometry.isperiodic

    natom = geometry.natom
    vr = geometry.vr

    # Box Vectors A = [a b c]
    A = Array(Diagonal(geometry.cell))
    T = typeof(A).parameters[begin]
    a, b, c = mapslices(x->[x], A, dims=[1])
    # shift vector
    vlmn = collect(Iterators.product(-1:1, -1:1, -1:1))
    vS = map(vlmn) do (rx, ry, rz)
        # 周期境界イメージのパターン, jの方だけイメージ分ループする
        A * T[rx; ry; rz]
    end

    # いったん初期化
    geometry.vflistexpand = Vector{Vector{Tuple{Int64, Tuple{Int64, Int64, Int64}}}}()
    vflistexpand = geometry.vflistexpand
    for i in 1:natom push!(vflistexpand, Vector{Tuple{Int64, Tuple{Int64, Int64, Int64}}}()) end
    R = geometry.R
    # 原子ペアの全組み合わせ
    pairs = combinations(1:natom, 2)
    # TODO parallelable
    for (i, j) in pairs
        @threads for (lmn, S) in collect(zip(vlmn, vS))
            vrⱼₛ = vr[:, j] + S
            r = norm(vr[:, i] - vrⱼₛ)
            if r < R
                push!(vflistexpand[i], (j, lmn))
            end
        end
    end
    return vflistexpand
end

"""
make_geometry

make example geometry
"""
function make_geometry(;isperiodic=false)
    # ジオメトリ作成
    nstep = 1000
    Δt = 1e-2
    cell = 10 * [1, 1, 1]
    atomnumbers = [1, 6]

    lj = LennardJones()
    ff = generalforcefunction(forcefunc(lj))
    forcemap = ForceMap(
        Dict(
            Pair(1, 1) => ff,
            Pair(1, 6) => ff,
            Pair(6, 6) => ff
        )
    )

    # 各原子数
    vnatom = [6, 12]
    natom = sum(vnatom)
    ntype = 2
    # [1, 1, 1, 1,  1,  1,  1,  1,  ... 2,  2,  2,  2,  2]
    vtype = vcat(map(x -> atomnumbers[x[1]]*ones(Int64, x[2]), enumerate(vnatom))...)
    vm = ones(natom)'
    vr = 10*init_vector((3, natom), DISTRIBUTE_UNIFORM)
    T = typeof(vr).parameters[1]
    A = Array(Diagonal{T}(cell))
    vr = modcoords(A, vr)
    # vv = zeros(3, natom)
    vv = 0*init_vector((3, natom), DISTRIBUTE_UNIFORM)
    # force
    R1 = 4
    R2 = 4
    R = 8
    L = 0
    vlist = Array{Array{Int64, 1}, 1}[]
    vlistexpand = Vector{Vector{Tuple{Int64, Tuple{Int64, Int64, Int64}}}}[]
    S = zeros(size(vr))
    for i in 1:natom push!(vlist, []) end
    vf0 = zeros(3, natom)
    vf = zeros(3, natom)
    return Geometry(
        nstep,
        1,
        Δt,
        cell,
        natom,
        atomnumbers,
        ntype,
        vtype,
        vm, vr, vv,
        R1, R2, R, L,
        vlist, vlistexpand, S,
        vf0, vf,
        repeat([isperiodic], 3),
        forcemap
    )
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

"""
step

calculate 1 step consider periodically
"""
function periodicalstep(geometry::Geometry)
    geometry.step += 1
    # periodically
    isloopx, isloopy, isloopz = geometry.isperiodic
    A = Array(Diagonal(geometry.cell))
    natom = geometry.natom
    Δt = geometry.Δt
    vm = geometry.vm
    vr = geometry.vr
    vv = geometry.vv
    vflistexpand = geometry.vflistexpand
    vf0 = geometry.vf0
    vf = geometry.vf
    # calculate step
    # cartesian coordinate
    vr .= vr .+ vv*Δt + (Δt^2 ./ 2vm .* vf)
    # printarray(vr)
    # fractional coordinate
    vrᵦ = fractionalcoords(A, vr)
    T = typeof(vrᵦ).parameters[1]

    # x, y, zで 0 ≤ p < 1 を満たさず， isloopの場合，座標を 0 ≤ p < 1に収まるようにする
    # rᵦᵢ = vrᵦ[:, i]
    cnt = 1
    IS = T[isloopx, isloopy, isloopz]
    # シフト行列 (3, natom) .* (3, 1)
    S = - floor.(vrᵦ) .* IS
    # printarray(S)
    # S = mapslices(vrᵦ, dims=1) do rᵦ
    #     # Matrix 
    #     floor.(vrᵦ)
    #     # periodicの領域外
    #     # 座標を丸める
    #     # -1.3 なら +2
    #     # 1.3　なら -1
    #     # 2.3 なら -2
    #     # floor
    #     return 
    # end
    # シフトして周期境界を適用する
    vrᵦ += S
    # printarray(vrᵦ)
    # printarray(vrᵦ)
    vr .= cartesiancoords(A, vrᵦ)
    # cartesian coordinate
    vf0 .= copy(vf)
    vf .= zeros(3, natom)

    # シフト行列を加算する
    geometry.S += S

    for (i, (vjlmn)) in enumerate(vflistexpand)
        @threads for (j, lmn) in vjlmn
            # shift vector
            s = Vector{T}([lmn...])
            # 原子間距離 シフト分の距離を戻す
            # r_12 = Float64[1.0, 1.0, 1.0]
            tᵢ = geometry.vtype[i]
            tⱼ = geometry.vtype[j]
            r_12 = 
            (vr[:, i] - cartesiancoords(A, geometry.S[:, i])) -
            (vr[:, j] - cartesiancoords(A, geometry.S[:, j]) + s)
            r = norm(r_12)

            nr_12 = r_12 / r
            force = geometry.forcemap(tᵢ, tⱼ)
            fvec = force(r_12...)
            vf[:, i] += fvec
            vf[:, j] -= fvec
        end
    end

    # calculate vv
    vv .= vv + (Δt^2 ./ 2vm .* (vf0 + vf))

    # Decrease Life of forcelistexpand
    geometry.Life += maximum(mapslices(norm, geometry.vr, dims=1)) * geometry.Δt
    if geometry.R2 < geometry.Life
        # @info "Recalculate forcelistexpand", geometry.step, geometry.Life
        forcelistexpand(geometry)
        # Reset Life and Shift-matrix
        geometry.S = zeros(size(geometry.vr))
        geometry.Life = 0
    end
end

"""
x̂ : fractional
x : cartesian
A : geometry
x = Ax̂
"""
function fractionalcoords(A::Matrix{T}, x::Matrix{T}) where T
    # Matrix{Float64}(Diagonal(cell))
    x̂ = inv(A) * x
    return x̂
end

function fractionalcoords(A::Matrix{T}, x::Vector{T}) where T
    # Matrix{Float64}(Diagonal(cell))
    x̂ = inv(A) * x
    return x̂
end

"""
x̂ : fractional
x : cartesian
A : geometry
x = Ax̂
"""
function cartesiancoords(A::Matrix{T}, x̂::Matrix{T}) where T
    x = A * x̂
    return x
end

function cartesiancoords(A::Matrix{T}, x̂::Vector{T}) where T
    x = A * x̂
    return x
end

"""
modcoords(A::Matrix{T}, x::Vector{T}) where T

atom position move in the cell
"""
function modcoords(A::Matrix{T}, x::Matrix{T}) where T
    x̂ = fractionalcoords(A, x)
    x̂ = mod.(x̂, 1)
    x = cartesiancoords(A, x̂)
end

function main(;logstep=20, isperiodic=true, seed=1)
    Random.seed!(seed)
    geom = make_geometry(isperiodic=isperiodic)
    forcelistexpand(geom)
    geoms = Geometry[copy(geom)]
    for i in 1:1000
        periodicalstep(geom)
        if i % logstep == 0
            push!(geoms, copy(geom))
        end
    end
    return geoms
end
