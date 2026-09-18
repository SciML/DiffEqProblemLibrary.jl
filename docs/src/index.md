# DiffEqProblemLibrary.jl

DiffEqProblemLibrary.jl is a component package in the DifferentialEquations
ecosystem. It contains premade problems for the differential equations solvers.
These can either be used as tests or as examples that show how to use the
solvers.

The problems are split into one sublibrary per problem class:

| Sublibrary | Problems |
| --- | --- |
| `ODEProblemLibrary` | ordinary differential equations |
| `SDEProblemLibrary` | stochastic differential equations |
| `DDEProblemLibrary` | delay differential equations |
| `DAEProblemLibrary` | differential-algebraic equations |
| `BVProblemLibrary` | boundary value problems |
| `JumpProblemLibrary` | jump and Gillespie problems |
| `NonlinearProblemLibrary` | nonlinear systems |

Each problem is exposed as a `prob_*` constructor, so a solver workbook can pull
in a canonical example instead of re-deriving one.

## Installation

To install DiffEqProblemLibrary.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("DiffEqProblemLibrary")
```

## Example

Each sublibrary is a dependency of the umbrella package, so one `using` per
sublibrary reaches its problems:

```@example problem-library
using DiffEqProblemLibrary
using ODEProblemLibrary: prob_ode_linear
using OrdinaryDiffEq: solve, Tsit5

sol = solve(prob_ode_linear, Tsit5())
sol.retcode
```

Users interested in the full solver API should see
[DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
and the [DifferentialEquations documentation](https://docs.sciml.ai/DiffEqDocs/stable/).
