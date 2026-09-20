using HomotopyContinuation
using Random

include(joinpath(@__DIR__, "AltsProblemLowerBound.jl"))

const LOWER_BOUND_SOLUTION_COUNT = 6 * 1442

"""
    certify_lower(; output_dir=joinpath(@__DIR__, "Certificates", "LowerBound"),
                  seed=20260918, max_precision=256, show_progress=true)

Find and certify all 8652 distinct solutions of `AltsProblem()` at an exact
real-rational parameter vector. Save the system and certified solution
intervals in `output_dir`, and return the certification result.
"""
function certify_lower(; output_dir::AbstractString=joinpath(@__DIR__, "Certificates", "LowerBound"),
                         seed::Integer=20260918, max_precision::Int=256,
                         show_progress::Bool=true)
    F = AltsProblem()
    expected = LOWER_BOUND_SOLUTION_COUNT

    println("Finding $expected lower-bound solutions by monodromy...")
    M = monodromy_solve(F; seed=UInt32(seed), target_solutions_count=expected,
                       min_solutions=expected, duplicate_check=:certified,
                       show_progress=show_progress)
    S = solutions(M)
    length(S) == expected ||
        error("Monodromy found $(length(S)) solutions; expected $expected.")

    P = parameters(M)
    rng = MersenneTwister(seed)
    nonzero_integer() = (rand(rng, Bool) ? 1 : -1) * rand(rng, 1:1000)
    Pnew = [BigInt(nonzero_integer()) // BigInt(nonzero_integer()) for _ in P]
    println("Tracking solutions to exact rational parameters...")
    T = solve(F, S; start_parameters=P, target_parameters=Pnew,
              show_progress=show_progress)
    tracked = [path.solution for path in path_results(T)]
    length(tracked) == expected ||
        error("Parameter tracking returned $(length(tracked)) paths; expected $expected.")
    all(!isnothing, tracked) || error("Some parameter paths have no endpoint.")

    println("Certifying $expected tracked solutions...")
    C = HomotopyContinuation.certify(F, tracked, Pnew;
                                     max_precision=max_precision,
                                     show_progress=show_progress)
    certified = ndistinct_certified(C)
    certified == expected ||
        error("Certified $certified distinct solutions; expected $expected.")

    println("Writing certificates to $output_dir...")
    mkpath(output_dir)
    open(joinpath(output_dir, "generic_system.txt"), "w") do io
        println(io, "system = AltsProblem()")
        println(io, "seed = ", seed)
        println(io, "expected = ", expected)
        println(io, "variables = ", variables(F))
        println(io, "parameters = ", HomotopyContinuation.parameters(F))
        println(io, "\nexact parameter values:")
        for (parameter, value) in zip(HomotopyContinuation.parameters(F), Pnew)
            println(io, parameter, " = ", value)
        end
        println(io, "\nequations:")
        for (i, equation) in enumerate(expressions(F))
            println(io, "F[", i, "] = ", equation)
        end
    end
    open(joinpath(output_dir, "certificates.txt"), "w") do io
        println(io, "system = AltsProblem()")
        println(io, "certified = ", certified)
        println(io, "\nCertified solution intervals:")
        for (i, certificate) in enumerate(distinct_certificates(C))
            println(io, "solution[", i, "] = ")
            show(io, certified_solution_interval(certificate))
            println(io)
        end
    end
    println("Saved $certified distinct certificates.")
    return C
end

"""
    certify_alt_mechanisms(; kwargs...)

Alias for `certify_lower` so the naming matches the paper's description.
"""
certify_alt_mechanisms(; kwargs...) = certify_lower(; kwargs...)
