"""
    normalize_boundary(B, H, s, u, epsilon, delta)

Construct S, L, M, polynomial K and factored guards for any registered regime.
Every route goes through finish_boundary_construction!, which preserves
L = M*H(S,delta^e) and makes K literally L(u,0) on the guarded chart.
Only the row operations are regime-specific; certification is shared.
"""
function normalize_boundary(B::BoundaryRegime, H, s, u, epsilon, delta)
    B.chart === nothing && error("$(B.id) has no recorded substitution chart")
    length(H) == length(s) == length(u) == length(B.omega) ||
        throw(ArgumentError("Boundary rows, original variables and local variables must match"))
    B.id === :B42 && return construct_B42(B, H, s, u, epsilon, delta)
    B.id === :B43 && return construct_B43(B, H, s, u, epsilon, delta)
    B.id === :B44 && return construct_B44(B, H, s, u, epsilon, delta)
    B.id in RECORDED_COMPATIBILITY_IDS &&
        return construct_recorded_boundary(B, H, s, u, epsilon, delta)
    B.id in GENERAL_COMPATIBILITY_IDS &&
        return construct_general_boundary(B, H, s, u, epsilon, delta)
    construct_scaled_boundary(B, H, s, u, epsilon, delta,
        (uq, dq) -> inv.(B.row_factors(uq, dq)); automatic_orders=true)
end

"""Construct the homotopy data automatically before normalizing a regime."""
function normalize_boundary(B::BoundaryRegime; seed::Integer=20260917)
    data = generate_homotopy(; seed=seed)
    normalize_boundary(B, data.H4, data.ss, data.u, data.epsilon, data.delta)
end
