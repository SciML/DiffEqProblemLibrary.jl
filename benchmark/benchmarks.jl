using Pkg
# Monorepo: the libraries live in lib/ — develop them so we benchmark this
# checkout's code, not the registered versions.
for lib in (
        "ODEProblemLibrary", "DAEProblemLibrary", "DDEProblemLibrary",
        "SDEProblemLibrary", "JumpProblemLibrary", "BVProblemLibrary",
        "NonlinearProblemLibrary",
    )
    Pkg.develop(path = joinpath(@__DIR__, "..", "lib", lib))
end

using ODEProblemLibrary, DAEProblemLibrary, NonlinearProblemLibrary
using BVProblemLibrary
using BenchmarkTools

const SUITE = BenchmarkGroup()

# =============================================================================
# Problem library loading (the package's job: constructing problems)
# =============================================================================

SUITE["ode"] = BenchmarkGroup()

SUITE["ode"]["vanderpol"] = @benchmarkable ODEProblemLibrary.prob_ode_vanderpol
SUITE["ode"]["lorenz"] = @benchmarkable ODEProblemLibrary.prob_ode_lorenz
SUITE["ode"]["rober"] = @benchmarkable ODEProblemLibrary.prob_ode_rober
SUITE["ode"]["brusselator_1d"] = @benchmarkable ODEProblemLibrary.prob_ode_brusselator_1d
SUITE["ode"]["pleiades"] = @benchmarkable ODEProblemLibrary.prob_ode_pleiades
SUITE["ode"]["hires"] = @benchmarkable ODEProblemLibrary.prob_ode_hires

SUITE["dae"] = BenchmarkGroup()

SUITE["dae"]["resrob"] = @benchmarkable DAEProblemLibrary.prob_dae_resrob

SUITE["nonlinear"] = BenchmarkGroup()

SUITE["nonlinear"]["testcase_lookup"] = @benchmarkable NonlinearProblemLibrary.nlprob_23_testcases["Generalized Rosenbrock function"]

SUITE["bvp"] = BenchmarkGroup()

SUITE["bvp"]["linear_1"] = @benchmarkable BVProblemLibrary.prob_bvp_linear_1
