const T = Float64
abstract type AbstractFilamentCache end
abstract type AbstractMagneticForce end
abstract type AbstractInextensibilityCache end
abstract type AbstractSolver end
abstract type AbstractSolverCache end
struct FerromagneticContinuous <: AbstractMagneticForce
    ω::T
    F::Vector{T}
end

mutable struct FilamentCache{
    MagneticForce <: AbstractMagneticForce,
    InextensibilityCache <: AbstractInextensibilityCache,
    SolverCache <: AbstractSolverCache
} <: AbstractFilamentCache
    N::Int
    μ::T
    Cm::T
    x::SubArray{T, 1, Vector{T}, Tuple{StepRange{Int, Int}}, true}
    y::SubArray{T, 1, Vector{T}, Tuple{StepRange{Int, Int}}, true}
    z::SubArray{T, 1, Vector{T}, Tuple{StepRange{Int, Int}}, true}
    A::Matrix{T}
    P::InextensibilityCache
    F::MagneticForce
    Sc::SolverCache
end
struct NoHydroProjectionCache <: AbstractInextensibilityCache
    J::Matrix{T}
    P::Matrix{T}
    J_JT::Matrix{T}
    J_JT_LDLT::LinearAlgebra.LDLt{T, SymTridiagonal{T}}
    P0::Matrix{T}

    function NoHydroProjectionCache(N::Int)
        new(zeros(N, 3 * (N + 1)),          # J
            zeros(3 * (N + 1), 3 * (N + 1)),    # P
            zeros(N, N),                 # J_JT
            LinearAlgebra.LDLt{T, SymTridiagonal{T}}(SymTridiagonal(zeros(N),
                zeros(N - 1))),
            zeros(N, 3 * (N + 1)))
    end
end
struct DiffEqSolverCache <: AbstractSolverCache
    S1::Vector{T}
    S2::Vector{T}

    DiffEqSolverCache(N::Integer) = new(zeros(T, 3 * (N + 1)), zeros(T, 3 * (N + 1)))
end
function FilamentCache(N = 20; Cm = 32, ω = 200, Solver = SolverDiffEq)
    InextensibilityCache = NoHydroProjectionCache
    SolverCache = DiffEqSolverCache
    tmp = zeros(3 * (N + 1))
    FilamentCache{FerromagneticContinuous, InextensibilityCache, SolverCache}(N, N + 1, Cm,
        view(tmp,
            1:3:(3 * (N + 1))),
        view(tmp,
            2:3:(3 * (N + 1))),
        view(tmp,
            3:3:(3 * (N + 1))),
        zeros(3 *
              (N + 1),
            3 *
            (N + 1)), # A
        InextensibilityCache(N), # P
        FerromagneticContinuous(ω,
            zeros(3 *
                  (N +
                   1))),
        SolverCache(N))
end
function stiffness_matrix!(f::AbstractFilamentCache)
    N, μ, A = f.N, f.μ, f.A
    A[:] = Matrix{Float64}(I, 3 * (N + 1), 3 * (N + 1))
    for i in 1:3
        A[i, i] = 1
        A[i, 3 + i] = -2
        A[i, 6 + i] = 1

        A[3 + i, i] = -2
        A[3 + i, 3 + i] = 5
        A[3 + i, 6 + i] = -4
        A[3 + i, 9 + i] = 1

        A[3 * (N - 1) + i, 3 * (N - 3) + i] = 1
        A[3 * (N - 1) + i, 3 * (N - 2) + i] = -4
        A[3 * (N - 1) + i, 3 * (N - 1) + i] = 5
        A[3 * (N - 1) + i, 3 * N + i] = -2

        A[3 * N + i, 3 * (N - 2) + i] = 1
        A[3 * N + i, 3 * (N - 1) + i] = -2
        A[3 * N + i, 3 * N + i] = 1

        for j in 2:(N - 2)
            A[3 * j + i, 3 * j + i] = 6
            A[3 * j + i, 3 * (j - 1) + i] = -4
            A[3 * j + i, 3 * (j + 1) + i] = -4
            A[3 * j + i, 3 * (j - 2) + i] = 1
            A[3 * j + i, 3 * (j + 2) + i] = 1
        end
    end
    rmul!(A, -μ^4)
    nothing
end
function update_separate_coordinates!(f::AbstractFilamentCache, r)
    N, x, y, z = f.N, f.x, f.y, f.z
    @inbounds for i in 1:length(x)
        x[i] = r[3 * i - 2]
        y[i] = r[3 * i - 1]
        z[i] = r[3 * i]
    end
    nothing
end

function update_united_coordinates!(f::AbstractFilamentCache, r)
    N, x, y, z = f.N, f.x, f.y, f.z
    @inbounds for i in 1:length(x)
        r[3 * i - 2] = x[i]
        r[3 * i - 1] = y[i]
        r[3 * i] = z[i]
    end
    nothing
end

function update_united_coordinates(f::AbstractFilamentCache)
    r = zeros(T, 3 * length(f.x))
    update_united_coordinates!(f, r)
    r
end

function initialize!(initial_conf_type::Symbol, f::AbstractFilamentCache)
    N, x, y, z = f.N, f.x, f.y, f.z
    if initial_conf_type == :StraightX
        x[:] = range(0, stop = 1, length = N + 1)
        y[:] = 0 .* x
        z[:] = 0 .* x
    else
        error("Unknown initial configuration requested.")
    end
    update_united_coordinates(f)
end

function magnetic_force!(::FerromagneticContinuous, f::AbstractFilamentCache, t)
    # TODO: generalize this for different magnetic fields as well
    N, μ, Cm, ω, F = f.N, f.μ, f.Cm, f.F.ω, f.F.F
    F[1] = -μ * Cm * cos(ω * t)
    F[2] = -μ * Cm * sin(ω * t)
    F[3 * (N + 1) - 2] = μ * Cm * cos(ω * t)
    F[3 * (N + 1) - 1] = μ * Cm * sin(ω * t)
    nothing
end

struct SolverDiffEq <: AbstractSolver end

function (f::FilamentCache)(dr, r, p, t)
    @views f.x, f.y, f.z = r[1:3:end], r[2:3:end], r[3:3:end]
    jacobian!(f)
    projection!(f)
    magnetic_force!(f.F, f, t)
    A, P, F, S1, S2 = f.A, f.P.P, f.F.F, f.Sc.S1, f.Sc.S2

    # implement dr = P * (A*r + F) in an optimized way to avoid temporaries
    mul!(S1, A, r)
    S1 .+= F
    mul!(S2, P, S1)
    copyto!(dr, S2)
    return dr
end

function (f::FilamentCache)(::Type{Val{:jac}}, J, r, p, t)
    mul!(J, f.P.P, f.A)
    nothing
end

function jacobian!(f::FilamentCache)
    N, x, y, z, J = f.N, f.x, f.y, f.z, f.P.J
    @inbounds for i in 1:N
        J[i, 3 * i - 2] = -2 * (x[i + 1] - x[i])
        J[i, 3 * i - 1] = -2 * (y[i + 1] - y[i])
        J[i, 3 * i] = -2 * (z[i + 1] - z[i])
        J[i, 3 * (i + 1) - 2] = 2 * (x[i + 1] - x[i])
        J[i, 3 * (i + 1) - 1] = 2 * (y[i + 1] - y[i])
        J[i, 3 * (i + 1)] = 2 * (z[i + 1] - z[i])
    end
    nothing
end

function projection!(f::FilamentCache)
    # implement P[:] = I - J'/(J*J')*J in an optimized way to avoid temporaries
    J, P, J_JT, J_JT_LDLT, P0 = f.P.J, f.P.P, f.P.J_JT, f.P.J_JT_LDLT, f.P.P0
    mul!(J_JT, J, J')
    LDLt_inplace!(J_JT_LDLT, J_JT)
    ldiv!(P0, J_JT_LDLT, J)
    mul!(P', P0, J)
    subtract_from_identity!(P)
    nothing
end

function subtract_from_identity!(A)
    rmul!(-1, A)
    @inbounds for i in 1:size(A, 1)
        A[i, i] += 1
    end
    nothing
end

function LDLt_inplace!(L::LinearAlgebra.LDLt{T, SymTridiagonal{T}},
        A::Matrix{T}) where {T <: Real}
    n = size(A, 1)
    dv, ev = L.data.dv, L.data.ev
    @inbounds for (i, d) in enumerate(diagind(A))
        dv[i] = A[d]
    end
    @inbounds for (i, d) in enumerate(diagind(A, -1))
        ev[i] = A[d]
    end
    @inbounds @simd for i in 1:(n - 1)
        ev[i] /= dv[i]
        dv[i + 1] -= abs2(ev[i]) * dv[i]
    end
    L
end

function filament_prob(::SolverDiffEq; N = 20, Cm = 32, ω = 200, time_end = 1.0)
    f = FilamentCache(N, Solver = SolverDiffEq, Cm = Cm, ω = ω)
    jac = (J, r, p, t) -> f(Val{:jac}, J, r, p, t)
    r0 = initialize!(:StraightX, f)
    stiffness_matrix!(f)
    prob = ODEProblem(ODEFunction(f, jac = jac), r0, (0.0, time_end))
end
"""
Filament PDE Discretization

Notebook: [Filament.ipynb](http://nbviewer.jupyter.org/github/JuliaDiffEq/DiffEqBenchmarks.jl/blob/master/StiffODE/Filament.ipynb)

In this problem is a real-world biological model from a paper entitled Magnetic dipole with a flexible tail as a self-propelling microdevice. It is a system of PDEs representing a Kirchhoff model of an elastic rod, where the equations of motion are given by the Rouse approximation with free boundary conditions.

"""
prob_ode_filament = filament_prob(SolverDiffEq())
