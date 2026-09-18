# Boundary 7, regime 1 (B71).

function S_B71(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    boundary5_chart(u9, 3 * u9 - 2 + delta^2 * u3, u1 * delta^(-4), u2 * delta,
        u4 * delta^4, u5 * delta^4, u6 * delta^4, u7 * delta^3, u8 * delta)
end

guard_B71(u) = boundary5_guard(u[9], 3 * u[9] - 2, u[1], u[2])

boundary_B71 = BoundaryRegime(
    "B_7, omega_7/4",
    [-1//1, -1//2, -3//2, -3//4, -5//4, -1//2, 1//4, -1//1, -1//4],
    4,
    [9],
    [1, 2, 3, 4, 5, 6, 7, 8],
    S_B71,
    _no_extra_factors,
    guard_B71
)
