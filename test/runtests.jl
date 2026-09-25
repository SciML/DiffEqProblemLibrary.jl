using Pkg
using SafeTestsets, Test
using SciMLTesting

const GROUP = current_group()
const LIB_DIR = joinpath(dirname(@__DIR__), "lib")

function activate_qa_env()
    qa_dir = joinpath(@__DIR__, "qa")
    return activate_group_env(qa_dir; parent = dirname(@__DIR__))
end

# The root umbrella package's own test: it re-exports each sublibrary, so the test
# is that the whole thing builds with no implicit or stale explicit imports. The old
# runtests.jl ran this body for GROUP=All, GROUP=Core, and GROUP=QA alike, so it is
# wired as the `qa` body, listed in `all`, and reached for GROUP=Core via an umbrella.
function qa_group()
    activate_qa_env()
    return @time @safetestset "ExplicitImports" include("qa/qa.jl")
end

@time begin
    # Route a requested monorepo sublibrary test group to its local package.
    base_group, test_group = detect_sublibrary_group(GROUP, LIB_DIR)

    if !isempty(base_group) && isdir(joinpath(LIB_DIR, base_group))
        Pkg.activate(joinpath(LIB_DIR, base_group))
        withenv("GROUP" => test_group) do
            Pkg.test(base_group, julia_args = ["--check-bounds=auto", "--compiled-modules=yes", "--depwarn=yes"], force_latest_compatible_version = false, allow_reresolve = true)
        end
    else
        # Root-package group dispatch. The previous runtests.jl ran the same
        # ExplicitImports QA body for GROUP=All, GROUP=Core, and GROUP=QA. `run_tests`
        # owns that routing: `qa` is the body, `all = ["QA"]` runs it under "All", and
        # the "Core" umbrella expands to "QA" so GROUP=Core runs it too.
        run_tests(;
            core = () -> nothing,
            qa = qa_group,
            all = ["QA"],
            umbrellas = Dict("Core" => ["QA"]),
        )
    end
end # @time
