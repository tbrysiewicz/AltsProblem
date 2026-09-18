# Boundary 3, regime 2 (B32).

function S_B32(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    denominator = delta^2 * u9
    [1 / denominator,
     (u2^2 + delta * u6) / denominator,
     u1 / denominator,
     (u2 + delta * u5) / denominator,
     (u1 * u2 + delta * u3) / denominator,
     (u2 * (u2 + delta * u5) + delta^2 * u7) / denominator,
     u2,
     (-u1 * (u2^2 + delta * u6) + delta * u4) / denominator,
     (- (u2 + delta * u5) * (u2^2 + delta * u6) + delta^2 * u8) / denominator]
end

guard_B32(u) = u[1] * u[2] * u[3] * u[9]

boundary_B32 = BoundaryRegime(
    "B_3 special",
    [-1//1, -1//1, -1//1, -1//1, -1//1, -1//1, 0//1, -1//1, -1//1],
    2,
    [1, 2],
    [3, 4, 5, 6, 7, 8, 9],
    S_B32,
    _no_extra_factors,
    guard_B32
)
