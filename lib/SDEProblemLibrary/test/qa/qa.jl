using SDEProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(SDEProblemLibrary; ambiguities = false, persistent_tasks = false)
end
