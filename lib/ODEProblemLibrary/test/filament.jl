using ODEProblemLibrary, LinearAlgebra, Test

# The right hand side is dr = P(r) * (A * r + F(t)), where P projects onto the
# tangent space of the inextensibility constraints, so every squared segment
# length must have zero rate of change for any state and any time.
function max_segment_rate(prob, N, u, t)
    du = similar(u)
    prob.f(du, u, prob.p, t)
    R, dR = reshape(u, 3, N + 1), reshape(du, 3, N + 1)
    rate = maximum(
        abs(2 * dot(R[:, i + 1] - R[:, i], dR[:, i + 1] - dR[:, i])) for i in 1:N
    )
    return rate / max(1.0, maximum(abs, du)), du
end

@testset "prob_ode_filament" begin
    prob = prob_ode_filament
    N = length(prob.u0) ÷ 3 - 1
    @test length(prob.u0) == 3 * (N + 1)

    for (u, t) in (
            (copy(prob.u0), 0.0),
            (copy(prob.u0), 0.013),
            (prob.u0 .+ 0.01 .* sinpi.(range(0, 2, length = length(prob.u0))), 0.013),
        )
        rate, du = max_segment_rate(prob, N, u, t)
        @test all(isfinite, du)
        @test rate < 1.0e-10
    end

    J = zeros(length(prob.u0), length(prob.u0))
    prob.f.jac(J, prob.u0, prob.p, 0.0)
    @test all(isfinite, J)
end
