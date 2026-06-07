using Pkg
using SafeTestsets, Test

const GROUP = get(ENV, "GROUP", "All")

@time begin
    # Detect sublibrary test groups.
    # GROUP can be a bare sublibrary name (Core test group) or
    # "{sublibrary}_{TEST_GROUP}" for any custom group (e.g., QA).
    # Sublibraries declare their groups in test/test_groups.toml (when they
    # differ from the default Core + QA).
    lib_dir = joinpath(dirname(@__DIR__), "lib")

    # Scan underscores right-to-left to find the longest matching sublibrary prefix.
    function _detect_sublibrary_group(group, lib_dir)
        isdir(joinpath(lib_dir, group)) && return (group, "Core")
        for i in length(group):-1:1
            if group[i] == '_' && isdir(joinpath(lib_dir, group[1:(i - 1)]))
                return (group[1:(i - 1)], group[(i + 1):end])
            end
        end
        return (group, "Core")
    end
    base_group, test_group = _detect_sublibrary_group(GROUP, lib_dir)

    if isdir(joinpath(lib_dir, base_group))
        Pkg.activate(joinpath(lib_dir, base_group))
        # On Julia < 1.11, the [sources] section in Project.toml is not supported.
        # Manually Pkg.develop local path dependencies so CI tests the PR branch code.
        # The transitive walk also develops each developed dependency's own [sources].
        if VERSION < v"1.11.0-DEV.0"
            developed = Set{String}()
            push!(developed, normpath(joinpath(lib_dir, base_group)))
            specs = Pkg.PackageSpec[]
            queue = [joinpath(lib_dir, base_group)]
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
    elseif GROUP == "All" || GROUP == "Core" || GROUP == "QA"
        # The root umbrella package's own test: it re-exports each sublibrary,
        # so the test is that the whole thing builds with no implicit or stale
        # explicit imports.
        @time @safetestset "ExplicitImports" begin
            using DiffEqProblemLibrary
            using ExplicitImports
            using Test
            @test check_no_implicit_imports(DiffEqProblemLibrary) === nothing
            @test check_no_stale_explicit_imports(DiffEqProblemLibrary) === nothing
        end
    end
end # @time
