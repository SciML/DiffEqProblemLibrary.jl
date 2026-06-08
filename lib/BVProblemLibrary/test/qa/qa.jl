using BVProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(BVProblemLibrary; ambiguities = false)
end
