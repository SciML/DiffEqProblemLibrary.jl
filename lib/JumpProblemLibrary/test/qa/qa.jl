using SciMLTesting, JumpProblemLibrary, Test

run_qa(
    JumpProblemLibrary; explicit_imports = true,
    aqua_kwargs = (; ambiguities = false, persistent_tasks = false)
)
