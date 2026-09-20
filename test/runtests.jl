using Test
using TOML
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
                    ScalingDimSector(charge=Dict("n"=>0.0), dims=[0.0, 1.0, 2.0]),
                    ScalingDimSector(charge=Dict("n"=>1.0), dims=[0.125, 1.125]),
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

        @testset "pretty-printing" begin
            entry = query_runs(; db_path=tmp, symmetry="Z2")[1]
            # Just check none of these error, and that the one-liners contain
            # the key numbers rather than a raw multi-line struct dump.
            @test occursin("χ=24", sprint(show, entry.params))
            @test occursin("n=1", sprint(show, entry.iterations[2].sectors[2]))
            @test occursin("c=0.4998", sprint(show, entry.iterations[2]))
            @test occursin("2 iters", sprint(show, entry))
        end

        @testset "final_iteration and find_sector" begin
            entry = query_runs(; db_path=tmp, symmetry="Z2")[1]
            @test final_iteration(entry).iteration == 20

            sec = find_sector(final_iteration(entry); n=1)
            @test sec !== nothing
            @test sec.dims == [0.125, 1.125]
            @test find_sector(final_iteration(entry); n=99) === nothing
        end

        @testset "find_sector is symmetry-agnostic" begin
            # Half-integer charges (e.g. O(2)/SU(2) spins) are stored as
            # their real value, not doubled.
            r = CFTResults(iteration=0, normalization=1.0, sectors=[
                ScalingDimSector(charge=Dict("j"=>0.5, "s"=>2.0), dims=[3.0]),
                ScalingDimSector(dims=[0.0]),  # trivial sector: empty charge
            ])
            @test find_sector(r; j=0.5, s=2).dims == [3.0]
            @test find_sector(r).dims == [0.0]  # no keywords = trivial charge
            @test occursin("j=1/2", sprint(show, r.sectors[1]))
        end

        @testset "CFTResults.dims (flattened across sectors)" begin
            r = query_runs(; db_path=tmp, symmetry="Z2")[1].iterations[2]
            # sectors are [n=0 -> [0.0,1.0,2.0]] and [n=1 -> [0.125,1.125]]
            @test r.dims == [0.0, 0.125, 1.0, 1.125, 2.0]
            @test :dims in propertynames(r)

            empty_r = CFTResults(iteration=0, normalization=0.0)
            @test empty_r.dims == Float64[]
        end

        @testset "trajectories and plateau_estimate" begin
            entry = query_runs(; db_path=tmp, symmetry="Z2")[1]
            @test central_charge_trajectory(entry) == [0.4998]

            traj = [1.0, 0.5, 0.501, 0.499, 0.5005, 0.5, 10.0]  # last point "blows up"
            est = plateau_estimate(traj; nwin=3, skipfrac=0.0)
            @test est.value ≈ 0.5 atol=0.01
            @test 7 ∉ est.range   # the blown-up last point must be excluded

            pc = plateau_central_charge(entry; nwin=1, skipfrac=0.0)
            @test pc.value ≈ 0.4998
        end

        @testset "export_csv" begin
            hits = query_runs(; db_path=tmp, symmetry="Z2")
            summary_out = joinpath(tmp, "summary.csv")
            export_csv(hits; out=summary_out)
            @test isfile(summary_out)
            summary_text = read(summary_out, String)
            @test occursin("BTRG", summary_text)
            @test length(readlines(summary_out)) == length(hits) + 1  # + header

            run_out = joinpath(tmp, "run.csv")
            export_csv(hits[1]; out=run_out)
            @test isfile(run_out)
            run_text = read(run_out, String)
            @test "n" in split(readlines(run_out)[1], ",") # charge column named after the real field
            @test occursin("0.125", run_text)
        end

        @testset "export_json" begin
            hits = query_runs(; db_path=tmp, symmetry="Z2")
            json_out = joinpath(tmp, "export.json")
            export_json(hits; out=json_out)
            @test isfile(json_out)
            text = read(json_out, String)

            # No real JSON parser in the stdlib to round-trip against, so
            # check structurally: balanced brackets and the expected content.
            @test count(==('{'), text) == count(==('}'), text)
            @test count(==('['), text) == count(==(']'), text)
            @test occursin("\"model\":\"phi4_real\"", text)
            @test occursin("\"n\":1", text)  # Z2 sector charge, from the real field name
            @test !occursin("NaN", text)     # NaN/missing must serialize as null, not raw NaN
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
                    infer_params = _ -> "x",   # simple form: just the model name
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

    @testset "displayed central charge is plateau-based, not literal-last" begin
        mktempdir() do tmp
            params = RunParameters(model="m", symmetry="s", algorithm="a",
                                    chi=1, K=1, mu0_sq=0.0, lambda=0.0)
            # Stable around 0.5 for 20 steps (a realistic trajectory length),
            # then blows up on the last two — the displayed value must
            # reflect the plateau (~0.5), not the literal last iteration.
            noise = 0.001 .* [1, -1, 0.5, -0.5, 0.3, -0.3, 0.2, -0.2, 0.1, -0.1,
                               0.4, -0.4, 0.15, -0.15, 0.25, -0.25, 0.05, -0.05, 0.35, -0.35]
            trajectory = vcat(0.5 .+ noise, [50.0, 100.0])
            iters = [
                CFTResults(iteration=i, normalization=1.0, central_charge=c)
                for (i, c) in enumerate(trajectory)
            ]
            entry = insert_run!(params, iters; db_path=tmp)

            @test occursin("c=0.50", sprint(show, entry))

            idx = TOML.parsefile(joinpath(tmp, "index.toml"))["entry"]
            @test idx[1]["central_charge"] ≈ 0.5 atol=0.01

            catalog = read(generate_catalog(; db_path=tmp), String)
            @test occursin("0.5000", catalog) || occursin("0.50", catalog)

            summary_out = joinpath(tmp, "summary.csv")
            export_csv([entry]; out=summary_out)
            summary_text = read(summary_out, String)
            csv_cc = parse(Float64, split(readlines(summary_out)[2], ",")[11])  # central_charge column
            @test csv_cc ≈ 0.5 atol=0.01
            @test !occursin("50.0", summary_text) && !occursin("100.0", summary_text)
        end
    end

    @testset "correct_field!" begin
        mktempdir() do tmp
            params = RunParameters(model="phi_complex", symmetry="U(1)", algorithm="LoopTNR",
                                    chi=32, K=10, mu0_sq=-0.5, lambda=0.5)
            iters = [CFTResults(iteration=0, normalization=1.0, central_charge=0.0)]
            e1 = insert_run!(params, iters; db_path=tmp)
            e2 = insert_run!(RunParameters(model="phi4_complex", symmetry="O(2)", algorithm="LoopTNR",
                                            chi=32, K=12, mu0_sq=-0.5, lambda=0.5),
                              iters; db_path=tmp)

            n = correct_field!(:model, "phi_complex", "phi4_complex"; db_path=tmp)
            @test n == 1

            fixed = load_entry(e1.id; db_path=tmp)
            @test fixed.params.model == "phi4_complex"
            @test fixed.params.symmetry == "U(1)"  # untouched fields preserved
            @test fixed.id == e1.id                 # same id, not a new entry

            untouched = load_entry(e2.id; db_path=tmp)
            @test untouched.params.model == "phi4_complex"  # was already correct, unaffected

            @test isempty(query_runs(; db_path=tmp, model="phi_complex"))
            @test length(query_runs(; db_path=tmp, model="phi4_complex")) == 2

            @test_throws ErrorException correct_field!(:chi, 1, 2; db_path=tmp)
        end
    end
end
