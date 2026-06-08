using JumpProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(JumpProblemLibrary; ambiguities = false, persistent_tasks = false)
end
