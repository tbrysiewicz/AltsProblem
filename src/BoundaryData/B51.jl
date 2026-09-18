# Boundary 5, regime 1 (B51).

# Boundary 5: three regimes

function S_B51(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    boundary5_chart(u8, u9, u1 * delta^(-10), u2 * delta^2,
        u3 * delta^6, u4 * delta^6, u5 * delta^6, u6 * delta^3, u7 * delta^3)
end

guard_B51(u) = boundary5_guard(u[8], u[9], u[1], u[2])

boundary_B51 = BoundaryRegime(
    "B_5, omega_5/3",
    [-5//3, -1//1, -8//3, -4//3, -7//3, -1//1, 1//3, -2//1, -2//3],
    6,
    [8, 9],
    [1, 2, 3, 4, 5, 6, 7],
    S_B51,
    _no_extra_factors,
    guard_B51
)
