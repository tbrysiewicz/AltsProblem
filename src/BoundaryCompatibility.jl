const RECORDED_COMPATIBILITY_IDS = (:B91, :B131, :B141, :B161, :B181)
const GENERAL_COMPATIBILITY_IDS = (:B151, :B171)

has_compatibility_constructor(B::BoundaryRegime) =
    !isempty(B.facial_indices) || B.id in RECORDED_COMPATIBILITY_IDS ||
    B.id in GENERAL_COMPATIBILITY_IDS

function evaluate_guard_factor(f, u, X, parameters=eltype(u)[], parameter_values=[])
    length(u) == length(X) || throw(ArgumentError("Guard coordinate count mismatch"))
    length(parameters) == length(parameter_values) ||
        throw(ArgumentError("Guard parameter count mismatch"))
    positions = var_index.(collect(u))
    parameter_positions = var_index.(collect(parameters))
    bound = Set(vcat(positions, parameter_positions))
    for a in exponents(f)
        any(a[j] != 0 for j in eachindex(a) if !(j in bound)) &&
            error("Guard has an unbound parameter; supply its exact parameter values")
    end
    # Embed the exact rational components in the interval arithmetic of X.
    values = [one(X[1])*real(p) + (one(X[1])*imag(p))*im
              for p in parameter_values]
    value = zero(X[1])
    for i in 1:length(f)
        a = exponent_vector(f, i)
        c = coeff(f, i)
        q = BigInt(numerator(c)) // BigInt(denominator(c))
        term = one(X[1]) * q
        for j in eachindex(u)
            term *= X[j]^Int(a[positions[j]])
        end
        for j in eachindex(values)
            term *= values[j]^Int(a[parameter_positions[j]])
        end
        value += term
    end
    return value
end

function construct_general_boundary(B::BoundaryRegime, H, s, u, epsilon, delta)
    out = BoundaryCompatibilityGeneral.construct_boundary_compatibility(
        B.id, H, s, u, epsilon, delta; chart=B.chart)
    B.e == out.e || error("$(B.id) ramification disagrees with its constructor")
    B.count == out.expected_count || error("$(B.id) count disagrees with its constructor")
    B.L, B.M, B.K = out.L, out.M, out.K
    B.orders = out.row_orders
    B.guard_factors = out.guard_factors
    B.guard = (X, parameters=eltype(u)[], parameter_values=[]) ->
        [evaluate_guard_factor(f, u, X, parameters, parameter_values)
         for f in out.guard_factors]
    return B.L, B.M, B.K, B.guard
end

function construct_recorded_boundary(B::BoundaryRegime, H, s, u, epsilon, delta)
    out = BoundaryCompatibilityExtras.construct_boundary_compatibility(
        B.id, H, s, u, epsilon, delta; chart=B.chart)
    B.e == out.e || error("$(B.id) ramification disagrees with its constructor")
    B.count == out.expected_count || error("$(B.id) count disagrees with its constructor")

    positions = Set(var_index.(u))
    for f in out.guard_factors, a in exponents(f)
        any(a[j] != 0 for j in eachindex(a) if !(j in positions)) &&
            error("$(B.id) guard depends on variables outside u")
    end

    B.L = out.L
    B.M = out.M
    B.K = out.K
    B.orders = out.row_orders
    B.guard_factors = out.guard_factors
    B.guard = X -> prod((evaluate_guard_factor(f, u, X)
                         for f in out.guard_factors); init=one(X[1]))
    return B.L, B.M, B.K, B.guard
end
