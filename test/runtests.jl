using Pkg
using SafeTestsets, Test
using SciMLTesting

const GROUP = current_group()
const LIB_DIR = joinpath(dirname(@__DIR__), "lib")

# The root umbrella package's own QA env activation. It is kept as an explicit
# helper (rather than `run_tests`'s `env =` group spec) so the activate/instantiate
# behavior — in particular NOT `Pkg.develop`-ing the root or transitively developing
# the env's `[sources]` on Julia < 1.11 — stays byte-for-byte identical to the
# previous hand-written runtests.jl.
function activate_qa_env()
    Pkg.activate(joinpath(@__DIR__, "qa"))
    return Pkg.instantiate()
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
    # Monorepo sublibrary routing. The root reads GROUP to pick a `lib/<sub>`
    # sublibrary, transitively develops its `[sources]` on Julia < 1.11, then
    # `Pkg.test`s it with the sub-group handed off via DIFFEQPROBLEMLIBRARY_TEST_GROUP.
    # This is kept as an explicit pre-step (rather than delegated to `run_tests`'s
    # built-in `lib_dir` path) so the sublibrary `Pkg.test` invocation — `julia_args`,
    # `force_latest_compatible_version = false`, `allow_reresolve = true` — stays
    # byte-for-byte identical to the previous runtests.jl. (`run_tests`'s sublibrary
    # path only passes `allow_reresolve`, which would silently drop `--depwarn=yes` /
    # `force_latest_compatible_version`.)
    base_group, test_group = detect_sublibrary_group(GROUP, LIB_DIR)

    if !isempty(base_group) && isdir(joinpath(LIB_DIR, base_group))
        Pkg.activate(joinpath(LIB_DIR, base_group))
        # On Julia < 1.11, the [sources] section in Project.toml is not supported.
        # Manually Pkg.develop local path dependencies so CI tests the PR branch code.
        # The transitive walk also develops each developed dependency's own [sources].
        if VERSION < v"1.11.0-DEV.0"
            developed = Set{String}()
            push!(developed, normpath(joinpath(LIB_DIR, base_group)))
            specs = Pkg.PackageSpec[]
            queue = [joinpath(LIB_DIR, base_group)]
            while !isempty(queue)
                pkg_dir = popfirst!(queue)
                toml_path = joinpath(pkg_dir, "Project.toml")
                isfile(toml_path) || continue
                toml = Pkg.TOML.parsefile(toml_path)
                if haskey(toml, "sources")
                    for (dep_name, source_spec) in toml["sources"]
                        if source_spec isa Dict && haskey(source_spec, "path")
                            dep_path = normpath(joinpath(pkg_dir, source_spec["path"]))
                            if isdir(dep_path) && !(dep_path in developed)
                                push!(developed, dep_path)
                                @info "Queuing local source dependency" dep_name dep_path
                                push!(specs, Pkg.PackageSpec(path = dep_path))
                                push!(queue, dep_path)
                            end
                        end
                    end
                end
            end
            isempty(specs) || Pkg.develop(specs)
        end
        withenv("DIFFEQPROBLEMLIBRARY_TEST_GROUP" => test_group) do
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
