using ODEProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(ODEProblemLibrary; ambiguities = false)
end
