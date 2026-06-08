using DDEProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(DDEProblemLibrary; ambiguities = false)
end
