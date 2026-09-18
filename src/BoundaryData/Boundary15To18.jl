# Boundaries 15--18 (B151, B161, B171, B181).

# Chart used for the exact certification of Boundaries 15--18.
function boundary15to18_chart(u1, u2, u3, u4, u5, u6, u7, u8, u9)
    u10 = (u3 - 1) * (u4 - 2 * u3 + 1) + u5
    u11 = -u3^2 + 2 * u3 + u9
    u12 = u3 * (u4 + 1 - u3) + u6
    u13 = -u10 * u3^2 + u8
    u14 = (u4 - 2 * u3 + 1) * u13 / u10 - u3^3 + u3^2 - u3 * u11 + u3 * u6 + u7
    [u1,
     u1 * u2^2 * u11,
     u1^2 * u2^2 * u10,
     u1 * u2 * u4,
     u1^2 * u2^3 * u3 * u10,
     u1 * u2^2 * u12,
     u2,
     u1^2 * u2^4 * u13,
     u1 * u2^3 * u14]
end

function scaled_boundary15to18_chart(coords, delta, powers)
    u1, u2, u3, u4, u5, u6, u7, u8, u9 = coords
    boundary15to18_chart(delta^powers[1] * u1, delta^powers[2] * u2, u3,
        delta^powers[3] * u4, delta^powers[4] * u5, delta^powers[5] * u6,
        delta^powers[6] * u7, delta^powers[7] * u8, delta^powers[8] * u9)
end

# Powers are for the eight scaled coordinates u1,u2,u4,...,u9; u3 has order zero.
function S_B151(u, delta)
    scaled_boundary15to18_chart(u, delta, (-8, 4, -2, 0, 0, -5, -7, -5))
end

function S_B161(u, delta)
    scaled_boundary15to18_chart(u, delta, (-2, 1, -1, 0, 0, -1, -2, -1))
end

function S_B171(u, delta)
    scaled_boundary15to18_chart(u, delta, (-4, 3, -3, -2, -2, -6, -7, -4))
end

function S_B181(u, delta)
    scaled_boundary15to18_chart(u, delta, (-3, 1, -2, 0, 0, 0, -2, 0))
end
