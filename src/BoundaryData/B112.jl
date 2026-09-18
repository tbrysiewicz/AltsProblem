# Boundary 11, regime 2 (B112).

# Compact local order: u8,u9 are free facial coordinates; u2 and u3 are
# the paper's u3 and u4 after eliminating the paper's u2.
function S_B112(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    eliminated = u2 * u3 / u8
    sigma = delta^2 * u1
    v8 = delta^7 * eliminated
    v4 = -2 * u8 * u9 + delta^4 * u2
    v2 = -2 + delta^3 * u3 + delta^4 * u4
    v5 = -u8 * u9 + delta^10 * u5
    v3 = -u8 * (1 - u8 * u9^2) - delta^4 * u8 * u2 * u9 + delta^10 * u6
    v7 = -1 + delta^3 * u3 + delta^4 * u4 - delta^7 * u8 * u9 * eliminated + delta^10 * u7
    sigma_chart(sigma, POWERS_B111_B112,
        [u8, v2, v3, v4, v5, one(u8), u9, v7, v8])
end

guard_B112(u) = u[8] * u[9] * u[1] * u[2] * u[3] * (1 - u[8] * u[9]^2)

boundary_B112 = BoundaryRegime(
    "B_11, regime 2",
    [7//5, 1//5, 8//5, 4//5, 1//1, 1//5, -3//5, 2//5, 3//10],
    10,
    [8, 9],
    [1, 2, 3, 4, 5, 6, 7],
    S_B112,
    _no_extra_factors,
    guard_B112
)
