# Boundary 4, regime 4 (B44).

function S_B44(u, delta)
    xi1, xi2, E, C, v, A, B, D, R = u
    eta = delta * E
    xi3 = xi2 + eta
    c = 2 * eta^2 + delta^3 * C
    d0 = xi2 * (xi2 + eta) * eta * c / (xi2 - eta)
    a = delta^5 * A
    b = delta^5 * B
    d = d0 + delta^5 * D
    e = -xi1 * c + d0 / xi2 + delta^5 * R
    sigma = delta^4 * v
    boundary4_chart(xi1, xi2, xi3, a, b, c, d, e, sigma)
end

guard_B44(u) = begin
    xi1, xi2, E, C, v, A, B, D, R = u
    v * E * xi1 * xi2 * (xi1 - xi2) * (xi1 - 2 * xi2) * (xi1 + 2 * xi2) * (xi1 + xi2)
end

boundary_B44 = BoundaryRegime(
    "B_4, regime 4",
    [-4//5, -4//5, -8//5, -4//5, -8//5, -4//5, 0//1, -8//5, -4//5],
    5,
    [1, 2],
    [3, 4, 5, 6, 7, 8, 9],
    S_B44,
    _no_extra_factors,
    guard_B44,
)

"""Construct L44_i = delta^-5 sigma^d_i H_i(S_B44(u, delta), delta^5)."""
function construct_B44(B::BoundaryRegime, H, s, u, epsilon, delta)
    d = (5, 5, 4, 5, 4, 4, 5, 4, 4)
    construct_scaled_boundary(B, H, s, u, epsilon, delta,
        (uQ, deltaQ) -> [deltaQ^(-5) * (deltaQ^4 * uQ[5])^d[i]
                         for i in eachindex(H)])
end
