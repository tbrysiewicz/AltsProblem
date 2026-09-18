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
        throw(ArgumentError("$regime requires the unresolved facial relation and further row cancellations; no constructor is supplied"))
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
