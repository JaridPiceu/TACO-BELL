using Test
using TACOBELL
using HDF5   # activates the ingest_jld2!/ingest_directory! extension

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

    @testset "ingest_directory! resilience" begin
        # ingest_directory! must not abort a hundred-file batch because one
        # file is unreadable — it should log and keep going. We don't have a
        # real TNRKit .jld2 fixture, so this exercises the error path with
        # deliberately-broken "jld2" files (plain text, not HDF5).
        mktempdir() do datadir
            mktempdir() do dbdir
                write(joinpath(datadir, "broken1.jld2"), "not an hdf5 file")
                write(joinpath(datadir, "broken2.jld2"), "also not an hdf5 file")
                mkpath(joinpath(datadir, "subdir"))
                write(joinpath(datadir, "subdir", "broken3.jld2"), "nope")

                summary = ingest_directory!(datadir;
                    infer_params = _ -> RunParameters(
                        model="x", symmetry="x", algorithm="x",
                        chi=0, K=0, mu0_sq=0.0, lambda=0.0,
                    ),
                    db_path = dbdir,
                )

                @test summary.total == 3
                @test summary.inserted == 0
                @test summary.failed == 3
                @test summary.skipped == 0
                @test isempty(query_runs(; db_path=dbdir))
            end
        end
    end
end
