# Boundary 4, regime 3 (B43).

function S_B43(u, delta)
    xi1, xi2, E, C, A, D, B, R, v = u
    sigma = delta * v
    xi3 = xi2 + sigma * E
    a = sigma^2 * A
    c = sigma * C
    d = sigma^2 * D
    b = xi1 * sigma^2 * A + sigma^3 * B
    e = -xi1 * sigma * C - xi1^2 * sigma^2 * A + sigma^2 * D / xi2 +
        sigma^3 * (R + xi1 * B)
    boundary4_chart(xi1, xi2, xi3, a, b, c, d, e, sigma)
end

guard_B43(u) = begin
    xi1, xi2, E, C, A, D, B, R, v = u
    v * E * C * xi1 * xi2 * (xi1 - xi2) * (xi1 - 2 * xi2) * (xi1 + 2 * xi2) * (xi1 + xi2)
end

boundary_B43 = BoundaryRegime(
    "B_4, regime 3",
    [-1//3, -1//3, -2//3, -1//3, -2//3, -1//3, 0//1, -2//3, -1//3],
    3,
    [1, 2],
    [3, 4, 5, 6, 7, 8, 9],
    S_B43,
    _no_extra_factors,
    guard_B43,
)

"""Construct L43_i = sigma^2 H_i(S_B43(u, delta), delta^3)."""
function construct_B43(B::BoundaryRegime, H, s, u, epsilon, delta)
    construct_scaled_boundary(B, H, s, u, epsilon, delta,
        (uQ, deltaQ) -> fill((deltaQ * uQ[9])^2, length(H)))
end
