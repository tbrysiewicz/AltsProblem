# Boundary 10, regime 1 (B101).

function S_B101(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    sigma = delta * u1
    u10 = -2 + delta * u3
    u11 = delta^3 * u2
    u12 = -2 * u8 * u9 + delta^2 * u4
    u13 = -u8 * u9 + delta^4 * u5
    u14 = -u8 * (1 - u8 * u9^2 + delta^2 * u4 * u9) + delta^4 * u6
    u15 = -1 + delta * u3 - delta^3 * u8 * u9 * u2 + delta^4 * u7
    sigma_chart(sigma, POWERS_B101, [u8, u10, u14, u12, u13, one(u8), u9, u15, u11])
end

guard_B101(u) = u[8] * u[9] * u[1] * u[2] * (1 - u[8] * u[9]^2)

boundary_B101 = BoundaryRegime(
    "B_10, regime 1",
    [5//4, 1//4, 3//2, 3//4, 1//1, 1//4, -1//2, 1//2, 1//2],
    4,
    [8, 9],
    [1, 2, 3, 4, 5, 6, 7],
    S_B101,
    _no_extra_factors,
    guard_B101
)
