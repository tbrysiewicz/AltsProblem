# Boundary 9, regime 1 (B91).

# The package constructs a full nine-variable compatibility system from this
# chart; no conversion to the archived eight-variable system is used.
function S_B91(u, delta)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = u
    u10, u11 = delta * u4, delta * u5
    u12 = u1 * u3 - u2 + delta * u6
    u13 = (u1 * u3 - u2) * u3 + delta * u7
    u14 = u2 - u1 * u3 + delta * u8
    u15 = -u1 + delta * u9
    [u11, 1 / u10, u10 * u11 * u12, u10 * u1, u10 * u11 * u13,
     u10 * u2, u3, u14, u15 / u11]
end
