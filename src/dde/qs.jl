"""
    prob_dde_qs

Delay differential equation model of Quorum Sensing (QS) of Pseudomonas putida IsoF in
continuous cultures.

# References

Buddrus-Schiemann et al. (2014). Analysis of N-Acylhomoserine Lactone Dynamics in Continuous
Cultures of Pseudomonas Putida IsoF By Use of ELISA and UHPLC/qTOF-MS-derived Measurements
and Mathematical Models, Analytical and Bioanalytical Chemistry.
"""
const prob_dde_qs =
  let S₀ = 1, τ = 2, D = 0.1, γₛ = 1.3e-12, Kₘ = 0.38, nₛ = 1.3, a = 0.66, αₐ = 2.3e-19,
    βₐ = 2.3e-18, C₁ = 70e-9, n₁ = 2.3, γₐ = 0.05, αC = 4e4, R = 5e-7, γ₃ = 0.08,
    Kₑ = 1.5e-4, αₗ = 1.1e-8, C₂ = 70e-9, n₂ = 2.5, γₗ = 0.005
    function hill_growth(x, K, n)
      if x > K
        return 1/(1 + (K / x)^n)
      end

      if K == 0
        return one(x)
      end

      y = (x / K)^n
      return y / (1 + y)
    end

    global function f_dde_qs!(du, u, h, p, t)
      if u[1] < 0
        # eq. 1 and 2 not defined for u[1] < 0
        du[1] = 0
        du[2] = 0
      else
        tmp = u[2] * hill_growth(u[1], Kₘ, nₛ)
        du[1] = D * (S₀ - u[1]) - γₛ * tmp
        du[2] = a * tmp - D * u[2]
      end

      if u[4] < 0
        # eq. 3 not defined for u[4] < 0
        du[3] = 0
        du[4] = max(0, αC * (R - u[4]) * u[3] - γ₃ * u[4])
      else
        du[3] = (αₐ + βₐ * hill_growth(u[4], C₁, n₁)) * u[2] -
          (γₐ + D) * u[3] - αC * (R - u[4]) * u[3] + γ₃ * u[4] -
          Kₑ * u[3] * u[5]
        du[4] = αC * (R - u[4]) * u[3] - γ₃ * u[4]
      end

      Cₜ = h(p, t - τ; idxs = 4)
      if Cₜ < 0
        # eq. 5 not defined for Cₜ < 0
        du[5] = 0
      else
        du[5] = αₗ * hill_growth(Cₜ, C₂, n₂) * u[2] - (γₗ + D) * u[5]
      end

      nothing
    end

    # history function for the fourth component
    function h_dde_qs(p, t; idxs = 4)
      idxs == 4 || error("history function is only implemented for the fourth component")
      t ≤ 0 || error("history function is only implemented for t ≤ 0")

      7.6e-8
    end

    DDEProblem(f_dde_qs!, [1, 8.4e8, 2.5e-9, 7.6e-8, 5e-15], h_dde_qs, (0.0, 45.0);
               constant_lags = [τ])
  end
