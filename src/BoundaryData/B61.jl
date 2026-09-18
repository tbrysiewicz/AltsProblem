# Boundary 6, regime 1 (B61).

# Boundaries 6--8

function S_B61(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    boundary5_chart(u8, u9, u1 * delta^(-6), u2 * delta,
        u3 * delta^4, u4 * delta^4, u5 * delta^4, u6 * delta^2, u7 * delta^2)
end

guard_B61(u) = boundary5_guard(u[8], u[9], u[1], u[2])

boundary_B61 = BoundaryRegime(
    "B_6, omega_6/4",
    [-3//2, -1//1, -5//2, -5//4, -9//4, -1//1, 1//4, -2//1, -3//4],
    4,
    [8, 9],
    [1, 2, 3, 4, 5, 6, 7],
    S_B61,
    _no_extra_factors,
    guard_B61
)
