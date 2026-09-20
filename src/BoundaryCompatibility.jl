module BoundaryCompatibilityExtras

import Oscar

export construct_B91, construct_B131, construct_B141, construct_B161, construct_B181,
       construct_boundary_compatibility

# Input: nine polynomial H4 rows in the manuscript's row order, all in one
# Oscar multivariate polynomial ring containing ss, u, epsilon and delta.
# Parameters remain unchanged by the substitution. Pass H4, not F4 alone.
# Output: L regular at delta=0, polynomial K, a rational matrix M, and a guard.
# L = M * H4(S(u,delta),delta^e), K = L(u,0).
# No numerical solve, certification, or expensive identity recheck is run.

function _generator_index(R, x)
    j = findfirst(==(x), Oscar.gens(R))
    isnothing(j) && throw(ArgumentError("Expected a generator of the input ring: $x"))
    return j
end

function _order(p, j)
    iszero(p) && throw(ArgumentError("An identically zero row needs a different construction"))
    return minimum(a[j] for a in Oscar.exponents(p))
end

function _valuation(f, j)
    return _order(Oscar.numerator(f), j) - _order(Oscar.denominator(f), j)
end

_at_zero(p, j) = Oscar.evaluate(p, [Oscar.gens(Oscar.parent(p))[j]],
                                [zero(Oscar.parent(p))])

function _boundary4_chart(x, y, z, a, b, c, d, r, sigma)
    [1 / sigma,
     (x^2 - 2*x*y + c) / sigma,
     y*z / sigma^2,
     (x + y + z + a) / sigma,
     x*y*z / sigma^2,
     (x*(x + z) + b) / sigma,
     x - y,
     (-x^2*y*z + d) / sigma^2,
     (-x^2*(x - y + z) + r) / sigma]
end

function _S_B91(u, delta)
    a, b = delta*u[4], delta*u[5]
    p = u[1]*u[3] - u[2]
    [b, 1/a, a*b*(p + delta*u[6]), a*u[1],
     a*b*(p*u[3] + delta*u[7]), a*u[2], u[3],
     -p + delta*u[8], (-u[1] + delta*u[9])/b]
end

function _S_B13_14(u, delta, m)
    sigma = delta^m*u[2]
    _boundary4_chart(delta*u[1], u[8], u[9],
        delta^m*u[3], delta^(m+1)*u[1]*u[4],
        delta^m*u[5], delta^m*u[6], delta^m*u[7], sigma)
end

function _S_B181(u, delta)
    T, W, tau = delta^(-3)*u[1], delta*u[2], u[3]
    h = delta^(-2)*u[4]
    a, b, c, d, e = u[5], u[6], u[7], delta^(-2)*u[8], u[9]
    Q = (tau - 1)*(h - 2*tau + 1) + a
    Z = -tau^2 + 2*tau + e
    N = tau*(h + 1 - tau) + b
    V = -Q*tau^2 + d
    R = (h - 2*tau + 1)*V/Q - tau^3 + tau^2 - tau*Z + tau*b + c
    [T, T*W^2*Z, T^2*W^2*Q, T*W*h, T^2*W^3*tau*Q,
     T*W^2*N, W, T^2*W^4*V, T*W^3*R]
end

function _spec(regime)
    regime == :B91 && return (1, _S_B91, 8)
    regime == :B131 && return (3, (u,d) -> _S_B13_14(u,d,3), 6)
    regime == :B141 && return (4, (u,d) -> _S_B13_14(u,d,4), 4)
    regime == :B161 && return (1, nothing, 2)
    regime == :B181 && return (2, _S_B181, 2)
    if regime in (:B151, :B171)
        throw(ArgumentError("$regime uses BoundaryCompatibilityGeneral; call normalize_boundary(boundary_regime(:$regime), ...)"))
    end
    throw(ArgumentError("Unsupported regime $regime"))
end

"""
    construct_boundary_compatibility(regime, H4, ss, u, epsilon, delta; chart=nothing)

Construct boundary systems directly from the supplied homotopy, retaining
all nine chart variables. `regime` is :B91, :B131, :B141, :B161, or :B181.
The optional `chart` callback must be the same substitution as the named chart.

For B131/B141, H4 must have the manuscript's unscaled row order. If input
rows have different constant scalings, undo those scalings before calling.

The returned `guard_factors` must all avoid zero on each certified interval.
`expected_count` is a manuscript target, not a certification result.
"""
function construct_boundary_compatibility(regime::Symbol, H4, ss, u, epsilon, delta;
                                          chart=nothing)
    length(H4) == length(ss) == length(u) == 9 ||
        throw(ArgumentError("Expected nine rows, nine s variables and nine u variables"))
    e, default_chart, expected = _spec(regime)
    regime == :B161 && isnothing(chart) &&
        throw(ArgumentError("B161 requires its recorded S_B161 chart"))
    R = Oscar.parent(H4[1])
    all(f -> Oscar.parent(f) === R, H4) ||
        throw(ArgumentError("All H4 rows must have the same polynomial parent"))
    Q = Oscar.fraction_field(R)
    ss, u = collect(ss), collect(u)
    indices = [_generator_index(R, x) for x in vcat(ss, u, [epsilon, delta])]
    length(unique(indices)) == 20 ||
        throw(ArgumentError("s, u, epsilon and delta must be distinct ring generators"))
    j = indices[end]
    dq, uq = Q(delta), Q.(u)
    S = Q.((isnothing(chart) ? default_chart : chart)(uq, dq))
    length(S) == 9 || throw(ArgumentError("Chart must return nine coordinates"))
    images = Q.(Oscar.gens(R))
    images[indices[1:9]] = S
    images[indices[end-1]] = dq^e
    substitute = Oscar.hom(R, Q, images)
    rows = [substitute(f) for f in H4]
    M = [i == k ? one(Q) : zero(Q) for i in 1:9, k in 1:9]

    if regime in (:B131, :B141)
        sigma, x = dq^e*uq[2], dq*uq[1]
        degrees = (5, 5, 4, 5, 4, 4, 5, 4, 4)
        for i in 1:9
            factor = sigma^degrees[i]
            rows[i] *= factor
            M[i,i] = factor
        end
        # Simultaneous operations on the Boundary 4 normalized rows.
        # Julia rows (2,4,7) and (1,4,7), not zero-based indices.
        rows[2] = rows[2] - 2*x*rows[4] + 3*x^2*rows[7]
        rows[1] = rows[1] - x^2*rows[4] + 2*x^3*rows[7]
        M[2,4], M[2,7] = -2*x*M[4,4], 3*x^2*M[7,7]
        M[1,4], M[1,7] = -x^2*M[4,4], 2*x^3*M[7,7]
    elseif regime == :B161
        T, W, tau, H, a, _, _, d, _ = uq
        Qscaled = (tau - 1)*(dq^(-1)*H - 2*tau + 1) + a
        initial_orders = (-6, -7, -6, -7, -6, -6, -7, -7, -6)
        for i in 1:9
            factor = dq^(-initial_orders[i])*Qscaled
            rows[i] *= factor
            M[i,i] = factor
        end
        # The row coefficients are the unscaled local variables, not chart values.
        rows[8] = (T*W^3*d*rows[8] + H*rows[2])/dq
        for k in 1:9
            M[8,k] = (T*W^3*d*M[8,k] + H*M[2,k])/dq
        end
    end

    # Guard each original leading coordinate, including its denominator.
    # This also covers the powers of u[2] used in the B13/B14 row scaling.
    factors = typeof(zero(R))[]
    function record!(p)
        iszero(p) && error("A required guard factor is identically zero")
        isone(p) || p in factors || push!(factors, p)
        return nothing
    end
    coordinate_orders = [_valuation(f, j) for f in S]
    for (f, a) in zip(S, coordinate_orders)
        leading_form = dq^(-a)*f
        record!(_at_zero(Oscar.numerator(leading_form), j))
        record!(_at_zero(Oscar.denominator(leading_form), j))
    end

    # First cancel exact rational expressions, then compute their valuations.
    # Never replace this step by the minimum weight of the unsummed monomials.
    row_orders = [_valuation(f, j) for f in rows]
    regime == :B161 && !all(iszero, row_orders) &&
        error("B161 rows are not regular after the recorded row operation")
    L = typeof(zero(R))[]
    K = typeof(zero(R))[]
    for i in 1:9
        f = dq^(-row_orders[i])*rows[i]
        n, d = Oscar.numerator(f), Oscar.denominator(f)
        record!(_at_zero(d, j))
        # Clear denominators in the full row, so K is literally L(delta=0).
        push!(L, n)
        push!(K, _at_zero(n, j))
        factor = Q(d)*dq^(-row_orders[i])
        for k in 1:9
            M[i,k] *= factor
        end
    end
    if regime == :B161
        T, W, H, d = u[1], u[2], u[4], u[8]
        monomials = (H*T^3*W^2, H*T^5*W^8*d, H*T^3*W^2,
                     H^3*T^5*W^7, H*T^3*W^2, H*T^3*W^2,
                     H^3*T^5*W^6, H*T^3*W^2, H*T^3*W^2)
        # These monomials divide K row by row. Divide the full L rationally;
        # the factors are nonzero on the guarded chart, including at delta=0.
        for x in (T, W, H, d)
            record!(x)
        end
        for i in 1:9
            K[i] = Oscar.divexact(K[i], monomials[i])
            for k in 1:9
                M[i,k] /= Q(monomials[i])
            end
        end
        L = [Q(L[i])/Q(monomials[i]) for i in 1:9]
    end
    guard = prod(factors; init=one(R))
    return (; L, M, K, guard, guard_factors=factors, S, e,
              row_orders, coordinate_orders, variables=u, expected_count=expected)
end

construct_B91(args...; kwargs...) = construct_boundary_compatibility(:B91, args...; kwargs...)
construct_B131(args...; kwargs...) = construct_boundary_compatibility(:B131, args...; kwargs...)
construct_B141(args...; kwargs...) = construct_boundary_compatibility(:B141, args...; kwargs...)
construct_B161(args...; kwargs...) = construct_boundary_compatibility(:B161, args...; kwargs...)
construct_B181(args...; kwargs...) = construct_boundary_compatibility(:B181, args...; kwargs...)

end # module

module BoundaryCompatibilityGeneral

import Oscar

export construct_B151, construct_B171, construct_boundary_compatibility

_at_zero(p, delta) = Oscar.evaluate(p, [delta], [zero(Oscar.parent(p))])

function _order(f, delta)
    iszero(f) && error("An identically zero expression requires another boundary chart")
    j = Oscar.var_index(delta)
    minimum(a[j] for a in Oscar.exponents(Oscar.numerator(f))) -
        minimum(a[j] for a in Oscar.exponents(Oscar.denominator(f)))
end

"""
    construct_boundary_compatibility(regime, H4, ss, u, epsilon, delta; chart)

Construct B151 or B171 for the package's unscaled, ordered homotopy
`H4 = (1-epsilon)*F4 + epsilon*G4`. Target parameters remain symbolic;
the coefficients of G4 are taken from H4. The existing nine-variable chart
is retained. No solve or certification is performed.

Returns rational rows L, a rational row matrix M, polynomial K, and guards,
with `L = M*H4(S,delta^e)` and literally `K = L(delta=0)` on the guarded
open set. The required parameter pivot is `ell[1,3] != 0`. Constants and
parameter specializations for which the recorded orders fail are rejected.
`row_orders` are the powers removed from the diagonal entries of M;
`initial_orders` and `row_divisions` record how they were obtained.
"""
function construct_boundary_compatibility(regime::Symbol, H4, ss, u, epsilon, delta;
                                          chart)
    regime in (:B151, :B171) || throw(ArgumentError("Unsupported regime $regime"))
    length(H4) == length(ss) == length(u) == 9 ||
        throw(ArgumentError("Expected nine rows, nine s variables and nine u variables"))
    A = Oscar.parent(H4[1])
    all(f -> Oscar.parent(f) === A, H4) ||
        throw(ArgumentError("All homotopy rows must have the same polynomial parent"))
    ss, u = collect(ss), collect(u)
    indices = Oscar.var_index.(vcat(ss, u, [epsilon, delta]))
    length(unique(indices)) == 20 ||
        throw(ArgumentError("s, u, epsilon and delta must be distinct generators"))
    Q = Oscar.fraction_field(A)
    dq = Q(delta)
    uq = Q.(u)
    T, W, tau, h, a, b, c, d, n = uq
    ramification = regime === :B151 ? 2 : 1
    expected_count = regime === :B151 ? 2 : 1
    initial_orders = regime === :B151 ?
        [-21, -24, -21, -23, -22, -21, -22, -23, -21] :
        [-13, -13, -13, -14, -13, -13, -13, -14, -13]
    qchart = regime === :B151 ?
        (tau-1)*(dq^(-2)*h-2*tau+1)+a :
        (tau-1)*(dq^(-3)*h-2*tau+1)+dq^(-2)*a

    # Only g_(2,2) contains s1^2*s9, with coefficient -1. Thus this
    # extracts ell[:,3] from F4 without assuming parameter names or values.
    F4 = [_at_zero(f, epsilon) for f in H4]
    ell3 = [-Oscar.coeff(f, ss, [2,0,0,0,0,0,0,0,1]) for f in F4]
    any(v -> any(a[Oscar.var_index(v)] != 0 for f in ell3
                 for a in Oscar.exponents(f)), vcat(ss, u, [epsilon, delta])) &&
        error("$regime expects target coefficients independent of the chart variables")
    pivot = ell3[1]
    iszero(pivot) && error("$regime requires the nonzero target pivot ell[1,3]")
    l = Q.(ell3)

    S = Q.(chart(uq, dq))
    images = Q.(Oscar.gens(A))
    images[indices[1:9]] = S
    images[indices[end-1]] = dq^ramification
    substitute = Oscar.hom(A, Q, images)
    scales = [qchart*dq^(-v) for v in initial_orders]
    rows = [scales[i]*substitute(H4[i]) for i in 1:9]
    M = [i == j ? scales[i] : zero(Q) for i in 1:9, j in 1:9]
    row_divisions = zeros(Int, 9)
    all(f -> _order(f, delta) == 0, rows) ||
        error("$regime input does not have the required initial row orders")

    function replace_row!(i, terms)
        rows[i] = sum(coef*rows[j] for (j, coef) in terms)/dq
        newrow = [sum(coef*M[j,k] for (j, coef) in terms)/dq for k in 1:9]
        M[i,:] = newrow
        row_divisions[i] += 1
        _order(rows[i], delta) == 0 ||
            error("$regime cancellation in row $i did not give order zero")
        return nothing
    end

    # Each operation is sequential. Multiplying by the pivot avoids
    # parameter denominators; the pivot is nevertheless a required guard.
    if regime === :B151
        replace_row!(2, ((2,l[1]), (1,T^2*W^5*d*(tau-1))))
        replace_row!(3, ((3,l[1]), (1,-l[3])))
        replace_row!(5, ((5,l[1]), (1,-T*W^3*n*(tau-1))))
        replace_row!(6, ((6,l[1]), (1,-l[6])))
        replace_row!(8, ((8,l[1]), (1,-h*T*W^2*(tau-1))))
        replace_row!(8, ((8,T*W^3*d), (2,h)))
        replace_row!(9, ((9,l[1]), (1,-l[9])))
        # The remaining leading relation is an identity of polynomials,
        # not merely an identity restricted to the leading solution curve.
        replace_row!(5, ((5,T*W^2*d), (2,n), (1,-T*W^2*d*l[5])))
    else
        replace_row!(3, ((3,l[1]), (1,-l[3])))
        replace_row!(5, ((5,l[1]), (1,-l[5])))
        replace_row!(6, ((6,l[1]), (1,-l[6])))
        replace_row!(8, ((8,l[1]), (1,-h*T*W^2*(tau-1))))
        replace_row!(9, ((9,l[1]), (1,-l[9])))
        replace_row!(2, ((2,l[1]), (1,T^2*W^5*d*(tau-1)-l[2])))
    end

    # These simple factors protect the row operations and all nine
    # prescribed leading s coordinates. a,b,c need not be nonzero.
    guard_factors = typeof(zero(A))[]
    function record!(p)
        iszero(p) && error("$regime has an identically zero guard factor")
        isone(p) || p in guard_factors || push!(guard_factors, p)
        return nothing
    end
    for p in vcat(u[[1,2,3,4,8,9]], [u[3]-1, pivot])
        record!(p)
    end
    coordinate_orders = [_order(f, delta) for f in S]
    for (f, v) in zip(S, coordinate_orders)
        leading = dq^(-v)*f
        record!(_at_zero(Oscar.numerator(leading), delta))
        record!(_at_zero(Oscar.denominator(leading), delta))
    end

    L = typeof(zero(Q))[]
    K = typeof(zero(A))[]
    for i in 1:9
        den0 = _at_zero(Oscar.denominator(rows[i]), delta)
        iszero(den0) && error("$regime row $i has a pole at delta=0")
        record!(den0)
        k = _at_zero(Oscar.numerator(rows[i]), delta)
        iszero(k) && error("$regime row $i needs another cancellation")
        scale = Q(den0)
        # Remove only factors already known to be units under the guard.
        # Apply the identical scaling to L and M, preserving both identities.
        for x in u[[1,2,3,4,8,9]]
            j = Oscar.var_index(x)
            power = minimum(ex[j] for ex in Oscar.exponents(k))
            if power > 0
                unit = x^power
                k = Oscar.divexact(k, unit)
                scale /= Q(unit)
            end
        end
        while iszero(Oscar.evaluate(k, [u[3]], [one(A)]))
            k = Oscar.divexact(k, u[3]-1)
            scale /= Q(u[3]-1)
        end
        push!(K, k)
        push!(L, scale*rows[i])
        for j in 1:9
            M[i,j] *= scale
        end
    end
    row_orders = initial_orders .+ row_divisions
    return (; L, M, K, S, guard_factors,
              guard=prod(guard_factors; init=one(A)),
              e=ramification, expected_count, variables=u,
              initial_orders, row_divisions, row_orders, coordinate_orders)
end

construct_B151(args...; kwargs...) = construct_boundary_compatibility(:B151, args...; kwargs...)
construct_B171(args...; kwargs...) = construct_boundary_compatibility(:B171, args...; kwargs...)

end # module


# B91 and B181 need only diagonal normalization and use the common route.
const RECORDED_COMPATIBILITY_IDS = (:B131, :B141, :B161)
const GENERAL_COMPATIBILITY_IDS = (:B151, :B171)

has_compatibility_constructor(B::BoundaryRegime) = B.chart !== nothing

function boundary_construction_method(B::BoundaryRegime)
    B.id in (:B131, :B141) && return :boundary4_row_cancellation
    B.id === :B161 && return :paired_row_cancellation
    B.id in GENERAL_COMPATIBILITY_IDS && return :parameter_row_cancellation
    B.id in (:B42, :B43, :B44) && return :prescribed_diagonal
    return :diagonal
end

_at_delta_zero(p, delta) = Oscar.evaluate(p, [delta], [zero(parent(p))])

function _record_guard!(factors, p)
    iszero(p) && error("A required boundary guard is identically zero")
    Oscar.total_degree(p) == 0 && return
    # Check factors individually instead of multiplying interval enclosures.
    for (f, _) in Oscar.factor(p)
        f in factors || push!(factors, f)
    end
end

function _record_leading_unit!(factors, f, delta)
    iszero(f) && error("A required boundary unit is identically zero")
    Q = parent(f)
    leading = f / Q(delta)^delta_order(f, delta)
    _record_guard!(factors, _at_delta_zero(numerator(leading), delta))
    _record_guard!(factors, _at_delta_zero(denominator(leading), delta))
end

"""Install the same exact construction contract for every boundary regime."""
function finish_boundary_construction!(B, out, H, s, u, epsilon, delta)
    A = parent(first(H)); Q = fraction_field(A)
    n = length(H)
    S = Q.(B.chart(Q.(u), Q(delta)))
    coordinate_orders = [delta_order(f, delta) for f in S]
    coordinate_orders == B.e .* B.omega ||
        error("$(B.id): chart orders $coordinate_orders disagree with e*omega=$(B.e .* B.omega)")
    factors = typeof(zero(A))[]
    for f in out.guard_factors
        _record_guard!(factors, f)
    end
    declared = Q(B.declared_guard(u))
    _record_guard!(factors, numerator(declared))
    _record_guard!(factors, denominator(declared))
    for f in S
        _record_leading_unit!(factors, f, delta)
    end

    M = zero_matrix(Q, n, n)
    for i in 1:n, j in 1:n
        M[i,j] = Q(out.M[i,j])
    end
    L = Q.(out.L)
    K = typeof(zero(A))[]
    for i in 1:n
        den0 = _at_delta_zero(denominator(L[i]), delta)
        iszero(den0) && error("$(B.id): row $i has a pole at delta=0")
        _record_guard!(factors, den0)
        limit = Q(_at_delta_zero(numerator(L[i]), delta)) / Q(den0)
        iszero(limit) && error("$(B.id): row $i needs another cancellation")
        # Taking numerator(limit) alone used to violate K=L(0).
        # Apply the identical denominator-clearing unit to L and M as well.
        unit = denominator(limit)
        _record_guard!(factors, unit)
        L[i] *= Q(unit)
        for j in 1:n
            M[i,j] *= Q(unit)
        end
        push!(K, numerator(limit))
    end
    # det(M) alone does not protect poles in individual matrix entries.
    for i in 1:n, j in 1:n
        iszero(M[i,j]) && continue
        den = denominator(M[i,j])
        leading = Q(den) / Q(delta)^delta_order(Q(den), delta)
        _record_guard!(factors, _at_delta_zero(numerator(leading), delta))
    end
    _record_leading_unit!(factors, det(M), delta)
    forbidden = var_index.(vcat(collect(s), [epsilon, delta]))
    for f in vcat(K, factors), a in exponents(f)
        any(j -> a[j] != 0, forbidden) &&
            error("$(B.id): a compatibility row or guard still contains s, epsilon or delta")
    end

    B.S, B.L, B.M, B.K = S, L, M, K
    B.orders = out.row_orders
    B.coordinate_orders = coordinate_orders
    B.guard_factors = factors
    B.guard = (X, parameters=eltype(u)[], parameter_values=[]) ->
        [evaluate_guard_factor(f, u, X, parameters, parameter_values) for f in factors]
    B.construction_metadata = (
        method=boundary_construction_method(B),
        initial_orders=get(out, :initial_orders, nothing),
        row_divisions=get(out, :row_divisions, nothing),
        coordinate_orders=coordinate_orders,
    )
    return L, M, K, B.guard
end

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
    return finish_boundary_construction!(B, out, H, s, u, epsilon, delta)
end

function construct_recorded_boundary(B::BoundaryRegime, H, s, u, epsilon, delta)
    out = BoundaryCompatibilityExtras.construct_boundary_compatibility(
        B.id, H, s, u, epsilon, delta; chart=B.chart)
    B.e == out.e || error("$(B.id) ramification disagrees with its constructor")
    B.count == out.expected_count || error("$(B.id) count disagrees with its constructor")

    return finish_boundary_construction!(B, out, H, s, u, epsilon, delta)
end
