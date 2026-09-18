# Boundary 8, regime 1 (B81).

function S_B81(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    u10 = u3 / delta
    u11 = -u3 * (u9 - 1) * u7 / delta + u8
    boundary5_chart(u9, u10, u1 * delta^(-3), u2 * delta, delta * u4, delta * u5, delta * u6,
        u11, u7)
end

function guard_B81(u)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    u1 * u2 * u9 * u3 * (u3 * u9 - u3) *
        (-u3 * u7 * u9 + u3 * u7 - u3 * u9^3 + u3 * u9^2) *
        (u7 - u9^2 + 2 * u9) * u8 * u7 * (u9 - 1)
end

boundary_B81 = BoundaryRegime(
    "B_8, omega_8/2",
    [-3//2, -1//2, -5//2, -3//2, -2//1, -1//1, 1//2, -3//2, -1//2],
    2,
    [9],
    [1, 2, 3, 4, 5, 6, 7, 8],
    S_B81,
    _no_extra_factors,
    guard_B81
)
