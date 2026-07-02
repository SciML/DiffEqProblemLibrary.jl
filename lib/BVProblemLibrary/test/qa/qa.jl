using SciMLTesting, BVProblemLibrary, Test

run_qa(BVProblemLibrary; explicit_imports = true, aqua_kwargs = (; ambiguities = false))
