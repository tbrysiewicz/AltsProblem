# Boundary 13, regime 1 (B131).

# Its compatibility constructor applies the manuscript's Boundary 4 row
# scaling and two row cancellations before taking leading coefficients.

function S_B131(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    sigma = delta^3 * u2
    boundary4_chart(delta * u1, u8, u9, delta^3 * u3, delta^4 * u1 * u4,
        delta^3 * u5, delta^3 * u6, delta^3 * u7, sigma)
end

guard_B131(u) = u[8] * u[9] * (u[9] - u[8]) * (u[9] + u[8]) * u[1] * u[2]
