# Boundary 2, regime 1 (B21).

function S_B21(u, delta)
    [delta * u[7], u[2], delta * u[7] * (-u[1] + delta * u[4]),
     delta * u[8], delta * u[7] * (-u[1] * u[3] + delta * u[5]),
     u[1], u[3], u[1] * (u[1] + u[2]) + delta * u[6], u[9] / u[7]]
end

factors_B21(u, delta) = fill(u[7], 9)
guard_B21(u) = u[1] * u[2] * u[3] * (u[1] + u[2]) * u[7] * u[8] * u[9]

boundary_B21 = BoundaryRegime(
    "B_2, regime 1",
    Rational{Int}.([1, 0, 1, 1, 1, 0, 0, 0, 0]),
    1,
    [1, 2, 3],
    [4, 5, 6, 7, 8, 9],
    S_B21,
    factors_B21,
    guard_B21
)
