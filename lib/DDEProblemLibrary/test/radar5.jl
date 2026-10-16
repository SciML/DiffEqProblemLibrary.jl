using DDEProblemLibrary, Test

@testset "RADAR5 Oregonator right-hand side" begin
    # F(1), F(2) of FCN in RADAR5-V2.1/OREGONATOR/dr-oregon.f
    k₁, k₂, k₃, k₄, f, A, B = 1.34, 1.6e9, 8.0e3, 4.0e7, 1.0, 6.0e-2, 6.0e-2
    y2_lag = 3.0e-6
    h(p, t; idxs = nothing) = idxs == 2 ? y2_lag : [0.0, y2_lag]
    u = [2.0e-8, 7.0e-6]
    expected = [
        k₁ * A * u[2] - k₂ * u[1] * y2_lag + k₃ * B * u[1] - 2 * k₄ * u[1]^2,
        -k₁ * A * u[2] - k₂ * u[1] * y2_lag + f * k₃ * B * u[1],
    ]
    prob = DDEProblemLibrary.prob_dde_RADAR5_oregonator
    du = zeros(2)
    prob.f(du, u, h, prob.p, 1.0)
    @test du ≈ expected rtol = 1.0e-12
end
