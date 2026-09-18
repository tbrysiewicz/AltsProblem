# Boundary 3, regime 1 (B31).

function S_B31(u, delta)
    [1 / (delta * u[9]), u[3] / (delta * u[9]), u[1] / (delta * u[9]),
     u[2] / (delta * u[9]), (u[1] * u[4] + delta * u[5]) / (delta * u[9]),
     (u[2] * u[4] + delta * u[6]) / (delta * u[9]), u[4],
     (-u[1] * u[3] + delta * u[7]) / (delta * u[9]),
     (-u[2] * u[3] + delta * u[8]) / (delta * u[9])]
end

guard_B31(u) = u[1] * u[2] * u[3] * u[4] * u[9] * (u[2] - u[4])

boundary_B31 = BoundaryRegime(
    "B_3 ordinary",
    Rational{Int}.([-1, -1, -1, -1, -1, -1, 0, -1, -1]),
    1,
    [1, 2, 3, 4],
    [5, 6, 7, 8, 9],
    S_B31,
    _no_extra_factors,
    guard_B31
)
