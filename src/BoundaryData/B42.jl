# Boundary 4, regime 2 (B42).

function S_B42(u, delta)
    xi1, xi2, xi3, A, C, D, B, R, v = u
    sigma = delta * v
    a = sigma * A
    c = sigma * C
    d = sigma * D
    b = xi1 * sigma * A + sigma^2 * B
    e = sigma * (-xi1^2 * A - xi1 * C + D / xi2) + sigma^2 * (R + xi1 * B)
    boundary4_chart(xi1, xi2, xi3, a, b, c, d, e, sigma)
end

guard_B42(u) = begin
    xi1, xi2, xi3, A, C, D, B, R, v = u
    v * xi1 * xi2 * xi3 * (xi3 - xi2) * (xi1 - xi2) * (xi1 - 2 * xi2) *
        (xi1 + xi2 + xi3) * (xi1 + xi3) * (xi1 - xi2 + xi3)
end

boundary_B42 = BoundaryRegime(
    "B_4, regime 2",
    [-1//2, -1//2, -1//1, -1//2, -1//1, -1//2, 0//1, -1//1, -1//2],
    2,
    [1, 2, 3],
    [4, 5, 6, 7, 8, 9],
    S_B42,
    _no_extra_factors,
    guard_B42,
)

"""Construct the exact B42 rows and their finite compatibility system."""
function construct_B42(B::BoundaryRegime, H, s, u, epsilon, delta)
    construct_scaled_boundary(B, H, s, u, epsilon, delta,
        (uQ, deltaQ) -> fill((deltaQ * uQ[9])^3, length(H)))
end
