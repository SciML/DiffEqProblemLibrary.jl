using DiffEqProblemLibrary
using Test
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems
using DiffEqProblemLibrary.DDEProblemLibrary: importddeproblems
using DiffEqProblemLibrary.DAEProblemLibrary: importdaeproblems
using DiffEqProblemLibrary.SDEProblemLibrary: importsdeproblems
using DiffEqProblemLibrary.JumpProblemLibrary: importjumpproblems

importodeproblems()
importddeproblems()
importjumpproblems()
importsdeproblems()
importdaeproblems()

# The test is simply that all of the examples build!
