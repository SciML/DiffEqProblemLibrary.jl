using NonlinearProblemLibrary
using Aqua

@testset "Aqua" begin
    Aqua.test_all(NonlinearProblemLibrary; ambiguities = false)
end
