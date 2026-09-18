# Boundary 12, regime 1 (B121).

function S_B121(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    sigma = delta * u1
    u10 = -2 + delta^2 * u2
    u11 = -u6 * u8 + delta^3 * u3
    u12 = -u6 * (u7 * u8 + u6 * u8^2 + delta * u1) + delta^3 * u4
    u13 = -1 - delta * u1 * u6 * u8 * u9 + delta^2 * u2 + delta^3 * u5
    sigma_chart(sigma, POWERS_B121,
        [u6, u10, u12, u7, u11, one(u6), u8, u13, u9])
end

guard_B121(u) = u[6] * u[7] * u[8] * u[9] * u[1] * (u[7] + u[6] * u[8])

boundary_B121 = BoundaryRegime(
    "B_12, regime 1",
    [2//1, 1//3, 2//1, 1//1, 4//3, 1//3, -1//1, 2//3, 0//1],
    3,
    [6, 7, 8, 9],
    [1, 2, 3, 4, 5],
    S_B121,
    _no_extra_factors,
    guard_B121
)
