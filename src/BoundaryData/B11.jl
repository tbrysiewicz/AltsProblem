# Boundary 1, regime 1 (B11).

function S_B11(u, delta)
    [delta * u[1], -2 * u[6] + delta * u[4],
     u[5]^2 - delta * (u[3] * u[5] + u[6] * u[1]) + delta^2 * u[8],
     -2 * u[5] + delta * u[3], -u[5] * u[6] + delta^2 * u[7],
     u[6], u[5] / (delta * u[1]),
     -u[6]^2 + delta * (u[6] * u[4] - u[2] * u[5]) + delta^2 * u[9],
     delta * u[2]]
end

guard_B11(u) = u[1] * u[2] * u[5] * u[6]

boundary_B11 = BoundaryRegime(
    "B_1, regime 1",
    Rational{Int}.([1//2, 0, 0, 0, 0, 0, -1//2, 0, 1//2]),
    2,
    [1, 2, 5, 6],
    [3, 4, 7, 8, 9],
    S_B11,
    _no_extra_factors,
    guard_B11
)
