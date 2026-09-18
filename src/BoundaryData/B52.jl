# Boundary 5, regime 2 (B52).

function S_B52(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    boundary5_chart(u9, 3 * u9 - 2 + delta^4 * u3, u1 * delta^(-10), u2 * delta^2,
        u4 * delta^10, u5 * delta^10, u6 * delta^10, u7 * delta^7, u8 * delta^3)
end

guard_B52(u) = boundary5_guard(u[9], 3 * u[9] - 2, u[1], u[2])

boundary_B52 = BoundaryRegime(
    "B_5, omega_5/5",
    [-1//1, -3//5, -8//5, -4//5, -7//5, -3//5, 1//5, -6//5, -2//5],
    10,
    [9],
    [1, 2, 3, 4, 5, 6, 7, 8],
    S_B52,
    _no_extra_factors,
    guard_B52
)
