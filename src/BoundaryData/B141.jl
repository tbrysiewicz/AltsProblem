# Boundary 14, regime 1 (B141).

# Boundary 14 substitution from the manuscript.
#
# Local order: u8 and u9 are the two free facial coordinates.

function S_B141(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    sigma = delta^4 * u2
    boundary4_chart(delta * u1, u8, u9, delta^4 * u3, delta^5 * u1 * u4,
        delta^4 * u5, delta^4 * u6, delta^4 * u7, sigma)
end

guard_B141(u) = guard_B131(u)
