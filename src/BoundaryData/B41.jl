# Boundary 4, regime 1 (B41).

# u1,u2,u3 are free facial coordinates; u4,...,u8 are first normal
# coefficients; u9 is the leading scale coefficient.

function S_B41(u, delta)
    [1 / (delta * u[9]),
     (u[1]^2 - 2 * u[2] * u[1] + delta * u[6]) / (delta * u[9]),
     u[2] * u[3] / (delta^2 * u[9]^2),
     (u[3] + u[1] + u[2] + delta * u[4]) / (delta * u[9]),
     u[1] * u[2] * u[3] / (delta^2 * u[9]^2),
     (u[1] * (u[3] + u[1]) + delta * u[5]) / (delta * u[9]),
     u[1] - u[2],
     (-u[2] * u[3] * u[1]^2 + delta * u[7]) / (delta^2 * u[9]^2),
     (-u[1]^2 * (u[3] + u[1] - u[2]) + delta * u[8]) / (delta * u[9])]
end

guard_B41(u) =
    u[1] * u[2] * u[3] * u[9] *
    (u[1] - 2 * u[2]) *
    (u[1] - u[2]) *
    (u[3] + u[1] + u[2]) *
    (u[3] + u[1]) *
    (u[3] + u[1] - u[2])

boundary_B41 = BoundaryRegime(
    "B_4, omega_4",
    Rational{Int}.([-1, -1, -2, -1, -2, -1, 0, -2, -1]),
    1,
    [1, 2, 3],
    [4, 5, 6, 7, 8, 9],
    S_B41,
    _no_extra_factors,
    guard_B41
)
