# Even if each polynomial is regular, the system need not be: leading terms can cancel between rows.

"""
    normalize_boundary(B, H, s, u, epsilon, delta)

Construct the exact row-normalized system for a boundary chart.  The current
operation records the diagonal row scaling made from the declared row factors
and the first nonzero power of `delta` in each row.  It returns `(L, M, K,
guard)`, where `L = M * H(S(u, delta))` and `K = L(u, 0)`.
"""
function normalize_boundary(B::BoundaryRegime, H, s, u, epsilon, delta)
    B.id === :B42 && return construct_B42(B, H, s, u, epsilon, delta)
    B.id === :B43 && return construct_B43(B, H, s, u, epsilon, delta)
    B.id === :B44 && return construct_B44(B, H, s, u, epsilon, delta)
    B.id in RECORDED_COMPATIBILITY_IDS &&
        return construct_recorded_boundary(B, H, s, u, epsilon, delta)
    B.id in GENERAL_COMPATIBILITY_IDS &&
        return construct_general_boundary(B, H, s, u, epsilon, delta)
    isempty(B.facial_indices) &&
        error("$(B.id) has no compatibility constructor")
    A = parent(H[1])
    Q = fraction_field(A)
    uQ = Q.(u)
    deltaQ = Q(delta)

    images = Q.(gens(A))
    images[var_index.(s)] = B.substitution(uQ, deltaQ)
    images[var_index(epsilon)] = deltaQ^B.e
    phi = hom(A, Q, images)
    substituted = phi.(H)

    extras = B.row_factors(uQ, deltaQ)
    rows = [substituted[i] / extras[i] for i in eachindex(substituted)]
    orders = [delta_order(row, delta) for row in rows]

    M = zero_matrix(Q, length(rows), length(rows))
    L = similar(rows)
    for i in eachindex(rows)
        scale = extras[i] * deltaQ^orders[i]
        M[i, i] = inv(scale)
        L[i] = rows[i] / deltaQ^orders[i]
    end

    K = map(L) do row
        numerator_at_zero = Oscar.evaluate(numerator(row), [delta], [zero(A)])
        denominator_at_zero = Oscar.evaluate(denominator(row), [delta], [zero(A)])
        @assert !iszero(denominator_at_zero)
        numerator(Q(numerator_at_zero) / Q(denominator_at_zero))
    end

    B.L = L
    B.M = M
    B.K = K
    B.orders = orders
    return L, M, K, B.guard
end

"""Construct the homotopy data automatically before normalizing a regime."""
function normalize_boundary(B::BoundaryRegime; seed::Integer=20260917)
    data = generate_homotopy(; seed=seed)
    normalize_boundary(B, data.H4, data.ss, data.u, data.epsilon, data.delta)
end
