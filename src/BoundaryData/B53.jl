# Boundary 5, regime 3 (B53).

function S_B53(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    boundary5_chart(u9, 3 * u9 - 2 + delta^3 * u3, u1 * delta^(-5), u2 * delta,
        u4 * delta^6, u5 * delta^6, u6 * delta^6, u7 * delta^5, u8 * delta^2)
end

guard_B53(u) = boundary5_guard(u[9], 3 * u[9] - 2, u[1], u[2])

boundary_B53 = BoundaryRegime(
    "B_5, omega_5/6",
    [-5//6, -1//2, -4//3, -2//3, -7//6, -1//2, 1//6, -1//1, -1//3],
    6,
    [9],
    [1, 2, 3, 4, 5, 6, 7, 8],
    S_B53,
    _no_extra_factors,
    guard_B53
)
