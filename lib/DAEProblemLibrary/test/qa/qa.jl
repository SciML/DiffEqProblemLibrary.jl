using SciMLTesting, DAEProblemLibrary, Test

run_qa(DAEProblemLibrary; explicit_imports = true, aqua_kwargs = (; ambiguities = false))
