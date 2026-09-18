# Boundary 11, regime 1 (B111).

# Compact local order: u7,u8,u9 are free facial coordinates; u2 is the
# paper's u3 after eliminating the paper's u2.
function S_B111(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    eliminated = (u8 + 2 * u7 * u9) * u2 / u7
    sigma = delta^2 * u1
    v8 = delta^3 * eliminated
    v2 = -2 + delta^3 * u2 + delta^4 * u3
    v5 = -u7 * u9 + delta^6 * u4
    v3 = -u7 * (1 + u8 * u9 + u7 * u9^2) + delta^6 * u5
    v7 = -1 + delta^3 * (u2 - u7 * u9 * eliminated) + delta^4 * u3 + delta^6 * u6
    sigma_chart(sigma, POWERS_B111_B112,
        [u7, v2, v3, u8, v5, one(u7), u9, v7, v8])
end

guard_B111(u) = u[7] * u[8] * u[9] * u[1] * u[2] *
    (u[8] + 2 * u[7] * u[9]) * (1 + u[8] * u[9] + u[7] * u[9]^2)

boundary_B111 = BoundaryRegime(
    "B_11, regime 1",
    [7//3, 1//3, 8//3, 4//3, 5//3, 1//3, -1//1, 2//3, -1//6],
    6,
    [7, 8, 9],
    [1, 2, 3, 4, 5, 6],
    S_B111,
    _no_extra_factors,
    guard_B111
)
