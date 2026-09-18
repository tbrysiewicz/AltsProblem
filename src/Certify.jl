module Certify

using HomotopyContinuation
using Arblib
using Oscar
using ..AltCertify: BoundaryRegime, generate_homotopy, finite_compatibility_system,
          hc_compatibility_system, has_compatibility_constructor,
          evaluate_guard_factor, RECORDED_COMPATIBILITY_IDS, GENERAL_COMPATIBILITY_IDS

export nearby_rational,
       nearby_rational_parameters,
       guard_excludes_zero,
       certify_regime,
       certify,
       reproduce_certification,
       CertificationRun

"""
Container for the complete numerical/exact certification pipeline.
"""
struct CertificationRun
    monodromy_result
    float_parameters
    soft_certification
    rational_parameters
    parameter_track_result
    tracked_solutions
    hard_certification
    guard_results
end

"""
    nearby_rational(x; tol=1e-12)

Return a nearby exact `Rational{BigInt}`.
"""
nearby_rational(x::Real; tol=1e-12) =
    rationalize(BigInt, x; tol=tol)

"""
    nearby_rational(z::Complex; tol=1e-12)

Return a nearby exact Gaussian rational.
"""
nearby_rational(z::Complex; tol=1e-12) =
    complex(nearby_rational(real(z); tol=tol),
            nearby_rational(imag(z); tol=tol))

nearby_rational_parameters(p; tol=1e-12) =
    [nearby_rational(z; tol=tol) for z in p]

# Current Arblib releases expose the underlying Arb predicate, but this
# wrapper tolerates the high-level spelling as well.
function _contains_zero(z)
    if isdefined(Arblib, :contains_zero)
        return getfield(Arblib, :contains_zero)(z)
    elseif isdefined(Arblib, :acb_contains_zero)
        return getfield(Arblib, :acb_contains_zero)(z) != 0
    elseif isdefined(Arblib, :arb_contains_zero)
        return getfield(Arblib, :arb_contains_zero)(z) != 0
    else
        error("Could not locate Arblib's contains-zero predicate.")
    end
end

"""
    guard_excludes_zero(guard, cert)

Evaluate the trusted guard on the certified solution box and prove that
its complex ball, or every ball in a returned list of guard factors,
does not contain zero.

`guard` accepts the certified solution interval (an `Arblib.AcbMatrix`)
and returns an `Acb`/`Arb` ball or a list of such balls.

Example:
    guard = X -> X[1] * X[2] * (X[1] + X[3])
"""
function guard_excludes_zero(guard::Function, cert)
    HomotopyContinuation.is_certified(cert) || return false
    X = HomotopyContinuation.certified_solution_interval(cert)
    y = guard(X)
    if y isa AbstractVector || y isa Tuple
        return all(z -> !_contains_zero(z), y)
    end
    return !_contains_zero(y)
end

function _require_distinct_certified(C, expected::Union{Nothing,Int}, label::String)
    n = HomotopyContinuation.ndistinct_certified(C)
    expected === nothing && return n
    n == expected || error("$label: certified $n distinct solutions; expected $expected.")
    return n
end

"""
    certify_regime(F, guard; expected=nothing, rational_tol=1e-12,
                   max_precision=256, monodromy_seed=nothing,
                   show_progress=true)

Run the certification pipeline for a parameterized square
`HomotopyContinuation.System` `F`:

1. Let monodromy choose a start pair and solve a Float64 fiber.
2. Soft-certify that numerical fiber.
3. Replace the parameter vector by a nearby exact Gaussian-rational vector.
4. Track all certified solutions to that nearby parameter vector.
5. Hard-certify at the exact rational parameter vector.
6. Evaluate the trusted guard on every certified solution interval and
   verify that zero is excluded.

The parameter homotopy uses a `ComplexF64` copy of the rational target,
but the final call to `certify` receives the exact rational vector.

`guard` is a function on the certified interval box, e.g.
    X -> X[1]*X[2]*X[5]*(X[1]+X[3])

If supplied, `guard_at_parameters(pstar)` binds the guard to the exact
parameter vector used by hard certification.

If `expected` is supplied, every stage is required to retain exactly that
many distinct certified solutions.
"""
function certify_regime(
    F,
    guard::Function;
    guard_at_parameters=nothing,
    expected::Union{Nothing,Int}=nothing,
    rational_tol::Real=1e-12,
    max_precision::Int=256,
    monodromy_seed::Union{Nothing,Integer}=nothing,
    show_progress::Bool=true,
)
    M = nothing
    sols = nothing
    for attempt in 1:10
    println("0. SEED THE MONODROMY SOLVE")
    #This is required for systems with positive dimensional components so that we 
    # stay on the one with regular solutions 

    Pinitial = randn(ComplexF64,63)
    Sinitial = solutions(HomotopyContinuation.solve(F;target_parameters = Pinitial))


    println("1. MONODROMY SOLVE")
    expected === nothing && error("certify_regime requires expected")
        M = HomotopyContinuation.monodromy_solve(
            F,Sinitial,Pinitial;
            target_solutions_count=expected,
            min_solutions=expected,
            unique_points_rtol=0.0000000000001,
            duplicate_check=:certified,
        )
        sols = HomotopyContinuation.solutions(M)
        length(sols) == expected && break
        attempt < 10 && println(
            "monodromy found $(length(sols)) of $(expected); retrying ($(attempt)/10)"
        )
    end
    pfloat = ComplexF64.(HomotopyContinuation.parameters(M))
    println("parameter count: ", length(pfloat))
    println("monodromy solutions: ", length(sols))
    expected === nothing || length(sols) == expected ||
        error("Monodromy found $(length(sols)) solutions; expected $expected.")

    println()
    println("2. SOFT CERTIFICATION AT FLOAT64 PARAMETERS")
    Csoft = HomotopyContinuation.certify(
        F, sols, HomotopyContinuation.parameters(M);
        max_precision=max_precision,
        show_progress=show_progress,
    )
    nsoft = _require_distinct_certified(Csoft, expected, "Soft certification")
    println("soft-certified distinct solutions: ", nsoft)

    println()
    println("3. NEARBY EXACT RATIONAL PARAMETERS")
    pstar = nearby_rational_parameters(pfloat; tol=rational_tol)
    maxmove = maximum(abs.(ComplexF64.(pstar) .- pfloat))
    println("maximum parameter displacement: ", maxmove)
    println("exact parameter vector constructed.")

    println()
    println("4. PARAMETER HOMOTOPY TO RATIONAL FIBER")
    T = HomotopyContinuation.solve(
        F, sols;
        start_parameters=HomotopyContinuation.parameters(M),
        target_parameters=pstar,
        show_progress=show_progress,
    )
    #Do not throw out points HC.jl thinks are singular
    tracked = [path.solution for path in HomotopyContinuation.path_results(T)]
    println("successfully tracked solutions: ", length(tracked))
    expected === nothing || length(tracked) == expected ||
        error("Parameter homotopy returned $(length(tracked)) solutions; expected $expected.")

    println()
    println("5. HARD CERTIFICATION AT EXACT RATIONAL PARAMETERS")
    Chard = HomotopyContinuation.certify(
        F, tracked, pstar;
        max_precision=max_precision,
        show_progress=show_progress,
    )
    nhard = _require_distinct_certified(Chard, expected, "Hard certification")
    println("hard-certified distinct solutions: ", nhard)

    println()
    println("6. GUARD CHECK ON CERTIFIED BOXES")
    hard_certs = HomotopyContinuation.distinct_certificates(Chard)
    exact_guard = isnothing(guard_at_parameters) ? guard : guard_at_parameters(pstar)
    guard_results = [guard_excludes_zero(exact_guard, c) for c in hard_certs]
    nguard = count(identity, guard_results)
    println("guard excludes zero on ", nguard, " / ", length(guard_results), " boxes")
    all(guard_results) ||
        error("Guard verification failed on $(length(guard_results)-nguard) certified boxes.")

    println()
    println("PASS: ", nhard, " distinct exact-rational certified solutions, all guarded.")

    return CertificationRun(
        M,
        pfloat,
        Csoft,
        pstar,
        T,
        tracked,
        Chard,
        guard_results,
    )
end

"""
    certify(B::BoundaryRegime; seed=20260917, output_root="Certificates", ...)

Construct and certify the compatibility system owned by `B`.  The exact
system, construction metadata, and interval certificates are written below
`output_root/B.id`.
"""
function certify(B::BoundaryRegime;
                 seed::Integer=20260917,
                 output_root::AbstractString="Certificates",
                 rational_tol::Real=1e-12,
                 max_precision::Int=256,
                 show_progress::Bool=true)
    B.chart === nothing && error("$(B.id) has no recorded substitution chart")
    has_compatibility_constructor(B) ||
        error("$(B.id) has no compatibility constructor")

    data = generate_homotopy(; seed=seed)
    K, _, orders = finite_compatibility_system(
        B, data.H4, data.ss, data.u, data.epsilon, data.delta)
    B.K = K
    B.orders = orders
    F = hc_compatibility_system(K, vec(data.ell);
        keep_parameters=true, seed=seed)
    guard_at_parameters = B.guard_factors === nothing ? nothing :
        p -> X -> [evaluate_guard_factor(f, data.u, X, vec(data.ell), p)
                   for f in B.guard_factors]
    run = certify_regime(F, B.guard;
        guard_at_parameters=guard_at_parameters,
        expected=B.count, rational_tol=rational_tol,
        max_precision=max_precision, monodromy_seed=seed,
        show_progress=show_progress)

    directory = joinpath(output_root, String(B.id))
    mkpath(directory)
    open(joinpath(directory, "generic_system.txt"), "w") do io
        println(io, "boundary = ", B.id)
        println(io, "seed = ", seed)
        println(io, "expected = ", B.count)
        println(io, "ramification = ", B.e)
        println(io, "delta_orders = ", orders)
        if B.id === :B42 || B.id === :B43 || B.id === :B44
            if B.id === :B42
                println(io, "normalization = L_i = (delta*u[9])^3 * H_i(S_B42(u,delta),delta^2)")
            elseif B.id === :B43
                println(io, "normalization = L_i = (delta*u[9])^2 * H_i(S_B43(u,delta),delta^3)")
            else
                println(io, "normalization = L_i = delta^-5 * (delta^4*u[5])^d_i * H_i(S_B44(u,delta),delta^5)")
                println(io, "d = (5,5,4,5,4,4,5,4,4)")
            end
            println(io, "chart = S_", B.id, ", defined in BoundaryData/", B.id, ".jl")
            println(io, "guard = ", factor(B.guard(data.u)))
        end
        if B.id in RECORDED_COMPATIBILITY_IDS || B.id in GENERAL_COMPATIBILITY_IDS
            constructor = B.id in GENERAL_COMPATIBILITY_IDS ?
                "BoundaryCompatibilityGeneral.construct_boundary_compatibility" :
                "BoundaryCompatibilityExtras.construct_boundary_compatibility"
            println(io, "construction = ", constructor)
            println(io, "row_orders = ", B.orders)
            println(io, "guard_factors =")
            for (i, f) in enumerate(B.guard_factors)
                println(io, "q[", i, "] = ", factor(f))
            end
            println(io, "S(u,delta) =")
            Q = fraction_field(parent(data.H4[1]))
            for (i, coordinate) in enumerate(B.chart(Q.(data.u), Q(data.delta)))
                println(io, "S[", i, "] = ", coordinate)
            end
            println(io, "M(u,delta) =")
            for i in 1:9, j in 1:9
                iszero(B.M[i, j]) || println(io, "M[", i, ",", j, "] = ", B.M[i, j])
            end
        end
        println(io, "\nlambda_star =")
        for (lambda, value) in zip(data.ell, run.rational_parameters)
            println(io, lambda, " = ", value)
        end
        println(io, "\nG4(s) =")
        for (i, equation) in enumerate(data.G4)
            println(io, "G4[", i, "] = ", equation)
        end
        println(io, "\nK(u; lambda) =")
        for (i, equation) in enumerate(K)
            println(io, "K[", i, "] = ", equation)
        end
    end
    open(joinpath(directory, "certificates.txt"), "w") do io
        println(io, "boundary = ", B.id)
        println(io, "certified = ", length(run.guard_results))
        println(io, "all_guards_exclude_zero = ", all(run.guard_results))
        println(io, "\nCertified solution intervals:")
        for (i, certificate) in enumerate(
            HomotopyContinuation.distinct_certificates(run.hard_certification))
            println(io, "solution[", i, "] = ")
            show(io, HomotopyContinuation.certified_solution_interval(certificate))
            println(io)
        end
    end
    return run
end

"""Re-run a saved certification, using the seed recorded in its system file."""
function reproduce_certification(B::BoundaryRegime;
                                 output_root::AbstractString="Certificates",
                                 kwargs...)
    directory = joinpath(output_root, String(B.id))
    system_file = joinpath(directory, "generic_system.txt")
    certificates_file = joinpath(directory, "certificates.txt")
    isfile(system_file) && isfile(certificates_file) ||
        error("no saved certification found in $directory")
    text = read(system_file, String)
    match_result = match(r"seed = (\d+)", text)
    seed = match_result === nothing ? 20260917 : parse(Int, match_result.captures[1])
    certify(B; output_root=output_root, seed=seed, kwargs...)
end

end # module
