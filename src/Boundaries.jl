_no_extra_factors(u, delta) = fill(one(u[1]), 9)

# ---------------------------------------------------------------------------
# Shared charts
# ---------------------------------------------------------------------------

function boundary4_chart(u1, u2, u3, eta1, eta2, eta3, eta4, eta5, sigma)
    [
        1/sigma,
        (u1^2 - 2*u1*u2 + eta3)/sigma,
        u2*u3/sigma^2,
        (u1 + u2 + u3 + eta1)/sigma,
        u1*u2*u3/sigma^2,
        (u1*(u1 + u3) + eta2)/sigma,
        u1 - u2,
        (-u1^2*u2*u3 + eta4)/sigma^2,
        (-u1^2*(u1 - u2 + u3) + eta5)/sigma,
    ]
end

function boundary5_chart(u8, u9, u10, u11, eta1, eta2, eta3, eta4, eta5)
    u1 = -u8^2 + 2*u8 + eta5
    u2 = (u8 - 1)*(u9 - 2*u8 + 1) + eta1
    u4 = u8*u2
    u5 = u8*(u9 + 1 - u8) + eta2
    u6 = -u2*u8^2 + eta4
    u7 = (u9 - 2*u8 + 1)*u6/u2 - u8^3 + u8^2 - u8*u1 + u8*eta2 + eta3
    [
        u10,
        u10*u11^2*u1,
        u10^2*u11^2*u2,
        u10*u11*u9,
        u10^2*u11^3*u4,
        u10*u11^2*u5,
        u11,
        u10^2*u11^4*u6,
        u10*u11^3*u7,
    ]
end

sigma_chart(sigma, powers, coordinates) =
    [sigma^powers[i]*coordinates[i] for i in 1:9]

function boundary5_guard(u8, u9, u10, u11)
    u1 = -u8^2 + 2*u8
    u2 = (u8 - 1)*(u9 - 2*u8 + 1)
    u4 = u8*u2
    u5 = u8*(u9 + 1 - u8)
    u6 = -u2*u8^2
    u7 = (u9 - 2*u8 + 1)*(-u8^2) - u8^3 + u8^2 - u8*u1
    u10*u11*u1*u2*u9*u4*u5*u6*u7
end
# ---------------------------------------------------------------------------
# Boundaries 10--12: sigma-toric charts
# ---------------------------------------------------------------------------

const POWERS_B101 = [5, 1, 6, 3, 4, 1, -2, 2, -1]
const POWERS_B111_B112 = [7, 1, 8, 4, 5, 1, -3, 2, -2]
const POWERS_B121 = [6, 1, 6, 3, 4, 1, -3, 2, 0]
const POWERS_B122 = [9, 1, 10, 5, 6, 1, -4, 2, -3]
