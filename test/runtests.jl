using Test
using TACOBELL

@testset "TACOBELL.jl" begin

    # Every test gets its own scratch database so runs don't leak between tests
    # or touch the real db/ directory.
    mktempdir() do tmp

        params = RunParameters(
            model="phi4_real", symmetry="Z2", algorithm="BTRG",
            chi=24, K=30, mu0_sq=-0.2, lambda=1.0,
        )

        iters = [
            CFTResults(
                iteration=0, normalization=0.0,
            ),
            CFTResults(
                iteration=20, normalization=0.693,
                central_charge=0.4998,
                sectors=[
                    ScalingDimSector(twice_j=0, s=0, dims=[0.0, 1.0, 2.0]),
                    ScalingDimSector(twice_j=2, s=1, dims=[0.125, 1.125]),
                ],
            ),
        ]

        @testset "insert and load" begin
            entry = insert_run!(params, iters; db_path=tmp)
            @test entry.params == params
            @test length(entry.iterations) == 2

            loaded = load_entry(entry.id; db_path=tmp)
            @test loaded.params == params
            @test length(loaded.iterations) == 2
            @test loaded.iterations[2].central_charge ≈ 0.4998
            @test length(loaded.iterations[2].sectors) == 2
        end

        @testset "duplicate protection" begin
            @test_throws ErrorException insert_run!(params, iters; db_path=tmp)
            entry2 = insert_run!(params, iters; db_path=tmp, allow_duplicate=true)
            @test entry2.id != ""
        end

        @testset "query_runs" begin
            hits = query_runs(; db_path=tmp, symmetry="Z2", chi=24)
            @test length(hits) == 2

            hits2 = query_runs(; db_path=tmp, chi=999)
            @test isempty(hits2)

            hits3 = query_runs(params; db_path=tmp)
            @test length(hits3) == 2
        end

        @testset "get_iteration" begin
            entry = query_runs(; db_path=tmp, symmetry="Z2")[1]
            r = get_iteration(entry, 20)
            @test r.central_charge ≈ 0.4998
            @test_throws ErrorException get_iteration(entry, 999)
        end

        @testset "find_closest" begin
            best = find_closest(-0.19, 1.0; db_path=tmp, symmetry="Z2", chi=24)
            @test best !== nothing
            @test best.params.mu0_sq == -0.2

            none = find_closest(0.0, 0.0; db_path=tmp, symmetry="does-not-exist")
            @test none === nothing
        end

        @testset "inspection" begin
            @test "BTRG" in list_algorithms(; db_path=tmp)
            @test "Z2" in list_symmetries(; db_path=tmp)
            # Just check it doesn't error.
            summarize_db(; db_path=tmp)
        end

        @testset "rebuild_index!" begin
            rm(joinpath(tmp, "index.toml"))
            rebuild_index!(; db_path=tmp)
            hits = query_runs(; db_path=tmp, symmetry="Z2")
            @test length(hits) == 2
        end

        @testset "generate_catalog" begin
            out = generate_catalog(; db_path=tmp)
            @test isfile(out)
            text = read(out, String)
            @test occursin("BTRG", text)
            @test occursin("2 run(s)", text)
        end
    end
end
