# Boundary 12, regime 2 (B122).

# Compact local order: u7,u8,u9 are free facial coordinates; u2 is the
# paper's u3 after eliminating the paper's u2.
function S_B122(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    eliminated = (u8 + 2 * u7 * u9) * u2 / u7
    sigma = delta * u1
    v8 = delta^2 * eliminated
    v2 = -2 + delta^2 * u2 + delta^3 * u3
    v5 = -u7 * u9 + delta^4 * u4
    v3 = -u7 * (1 + u8 * u9 + u7 * u9^2) + delta^4 * u5
    v7 = -1 + delta^2 * (u2 - u7 * u9 * eliminated) + delta^3 * u3 + delta^4 * u6
    sigma_chart(sigma, POWERS_B122,
        [u7, v2, v3, u8, v5, one(u7), u9, v7, v8])
end

guard_B122(u) = u[7] * u[8] * u[9] * u[1] * u[2] *
    (u[8] + 2 * u[7] * u[9]) * (1 + u[8] * u[9] + u[7] * u[9]^2)

boundary_B122 = BoundaryRegime(
    "B_12, regime 2",
    [9//4, 1//4, 5//2, 5//4, 3//2, 1//4, -1//1, 1//2, -1//4],
    4,
    [7, 8, 9],
    [1, 2, 3, 4, 5, 6],
    S_B122,
    _no_extra_factors,
    guard_B122
)
