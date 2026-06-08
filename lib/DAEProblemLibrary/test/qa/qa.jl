using DAEProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(DAEProblemLibrary; ambiguities = false)
end
