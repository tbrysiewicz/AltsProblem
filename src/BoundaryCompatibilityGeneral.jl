module BoundaryCompatibilityGeneral

import Oscar

export construct_B151, construct_B171, construct_boundary_compatibility

_at_zero(p, delta) = Oscar.evaluate(p, [delta], [zero(Oscar.parent(p))])

function _order(f, delta)
    iszero(f) && error("An identically zero expression requires another boundary chart")
    j = Oscar.var_index(delta)
    minimum(a[j] for a in Oscar.exponents(Oscar.numerator(f))) -
        minimum(a[j] for a in Oscar.exponents(Oscar.denominator(f)))
end

"""
    construct_boundary_compatibility(regime, H4, ss, u, epsilon, delta; chart)

Construct B151 or B171 for the package's unscaled, ordered homotopy
`H4 = (1-epsilon)*F4 + epsilon*G4`. Target parameters remain symbolic;
the coefficients of G4 are taken from H4. The existing nine-variable chart
is retained. No solve or certification is performed.

Returns rational rows L, a rational row matrix M, polynomial K, and guards,
with `L = M*H4(S,delta^e)` and literally `K = L(delta=0)` on the guarded
open set. The required parameter pivot is `ell[1,3] != 0`. Constants and
parameter specializations for which the recorded orders fail are rejected.
`row_orders` are the powers removed from the diagonal entries of M;
`initial_orders` and `row_divisions` record how they were obtained.
"""
function construct_boundary_compatibility(regime::Symbol, H4, ss, u, epsilon, delta;
                                          chart)
    regime in (:B151, :B171) || throw(ArgumentError("Unsupported regime $regime"))
    length(H4) == length(ss) == length(u) == 9 ||
        throw(ArgumentError("Expected nine rows, nine s variables and nine u variables"))
    A = Oscar.parent(H4[1])
    all(f -> Oscar.parent(f) === A, H4) ||
        throw(ArgumentError("All homotopy rows must have the same polynomial parent"))
    ss, u = collect(ss), collect(u)
    indices = Oscar.var_index.(vcat(ss, u, [epsilon, delta]))
    length(unique(indices)) == 20 ||
        throw(ArgumentError("s, u, epsilon and delta must be distinct generators"))
    Q = Oscar.fraction_field(A)
    dq = Q(delta)
    uq = Q.(u)
    T, W, tau, h, a, b, c, d, n = uq
    ramification = regime === :B151 ? 2 : 1
    expected_count = regime === :B151 ? 2 : 1
    initial_orders = regime === :B151 ?
        [-21, -24, -21, -23, -22, -21, -22, -23, -21] :
        [-13, -13, -13, -14, -13, -13, -13, -14, -13]
    qchart = regime === :B151 ?
        (tau-1)*(dq^(-2)*h-2*tau+1)+a :
        (tau-1)*(dq^(-3)*h-2*tau+1)+dq^(-2)*a

    # Only g_(2,2) contains s1^2*s9, with coefficient -1. Thus this
    # extracts ell[:,3] from F4 without assuming parameter names or values.
    F4 = [_at_zero(f, epsilon) for f in H4]
    ell3 = [-Oscar.coeff(f, ss, [2,0,0,0,0,0,0,0,1]) for f in F4]
    any(v -> any(a[Oscar.var_index(v)] != 0 for f in ell3
                 for a in Oscar.exponents(f)), vcat(ss, u, [epsilon, delta])) &&
        error("$regime expects target coefficients independent of the chart variables")
    pivot = ell3[1]
    iszero(pivot) && error("$regime requires the nonzero target pivot ell[1,3]")
    l = Q.(ell3)

    S = Q.(chart(uq, dq))
    images = Q.(Oscar.gens(A))
    images[indices[1:9]] = S
    images[indices[end-1]] = dq^ramification
    substitute = Oscar.hom(A, Q, images)
    scales = [qchart*dq^(-v) for v in initial_orders]
    rows = [scales[i]*substitute(H4[i]) for i in 1:9]
    M = [i == j ? scales[i] : zero(Q) for i in 1:9, j in 1:9]
    row_divisions = zeros(Int, 9)
    all(f -> _order(f, delta) == 0, rows) ||
        error("$regime input does not have the required initial row orders")

    function replace_row!(i, terms)
        rows[i] = sum(coef*rows[j] for (j, coef) in terms)/dq
        newrow = [sum(coef*M[j,k] for (j, coef) in terms)/dq for k in 1:9]
        M[i,:] = newrow
        row_divisions[i] += 1
        _order(rows[i], delta) == 0 ||
            error("$regime cancellation in row $i did not give order zero")
        return nothing
    end

    # Each operation is sequential. Multiplying by the pivot avoids
    # parameter denominators; the pivot is nevertheless a required guard.
    if regime === :B151
        replace_row!(2, ((2,l[1]), (1,T^2*W^5*d*(tau-1))))
        replace_row!(3, ((3,l[1]), (1,-l[3])))
        replace_row!(5, ((5,l[1]), (1,-T*W^3*n*(tau-1))))
        replace_row!(6, ((6,l[1]), (1,-l[6])))
        replace_row!(8, ((8,l[1]), (1,-h*T*W^2*(tau-1))))
        replace_row!(8, ((8,T*W^3*d), (2,h)))
        replace_row!(9, ((9,l[1]), (1,-l[9])))
        # The remaining leading relation is an identity of polynomials,
        # not merely an identity restricted to the leading solution curve.
        replace_row!(5, ((5,T*W^2*d), (2,n), (1,-T*W^2*d*l[5])))
    else
        replace_row!(3, ((3,l[1]), (1,-l[3])))
        replace_row!(5, ((5,l[1]), (1,-l[5])))
        replace_row!(6, ((6,l[1]), (1,-l[6])))
        replace_row!(8, ((8,l[1]), (1,-h*T*W^2*(tau-1))))
        replace_row!(9, ((9,l[1]), (1,-l[9])))
        replace_row!(2, ((2,l[1]), (1,T^2*W^5*d*(tau-1)-l[2])))
    end

    # These simple factors protect the row operations and all nine
    # prescribed leading s coordinates. a,b,c need not be nonzero.
    guard_factors = typeof(zero(A))[]
    function record!(p)
        iszero(p) && error("$regime has an identically zero guard factor")
        isone(p) || p in guard_factors || push!(guard_factors, p)
        return nothing
    end
    for p in vcat(u[[1,2,3,4,8,9]], [u[3]-1, pivot])
        record!(p)
    end
    coordinate_orders = [_order(f, delta) for f in S]
    for (f, v) in zip(S, coordinate_orders)
        leading = dq^(-v)*f
        record!(_at_zero(Oscar.numerator(leading), delta))
        record!(_at_zero(Oscar.denominator(leading), delta))
    end

    L = typeof(zero(Q))[]
    K = typeof(zero(A))[]
    for i in 1:9
        den0 = _at_zero(Oscar.denominator(rows[i]), delta)
        iszero(den0) && error("$regime row $i has a pole at delta=0")
        record!(den0)
        k = _at_zero(Oscar.numerator(rows[i]), delta)
        iszero(k) && error("$regime row $i needs another cancellation")
        scale = Q(den0)
        # Remove only factors already known to be units under the guard.
        # Apply the identical scaling to L and M, preserving both identities.
        for x in u[[1,2,3,4,8,9]]
            j = Oscar.var_index(x)
            power = minimum(ex[j] for ex in Oscar.exponents(k))
            if power > 0
                unit = x^power
                k = Oscar.divexact(k, unit)
                scale /= Q(unit)
            end
        end
        while iszero(Oscar.evaluate(k, [u[3]], [one(A)]))
            k = Oscar.divexact(k, u[3]-1)
            scale /= Q(u[3]-1)
        end
        push!(K, k)
        push!(L, scale*rows[i])
        for j in 1:9
            M[i,j] *= scale
        end
    end
    row_orders = initial_orders .+ row_divisions
    return (; L, M, K, S, guard_factors,
              guard=prod(guard_factors; init=one(A)),
              e=ramification, expected_count, variables=u,
              initial_orders, row_divisions, row_orders, coordinate_orders)
end

construct_B151(args...; kwargs...) = construct_boundary_compatibility(:B151, args...; kwargs...)
construct_B171(args...; kwargs...) = construct_boundary_compatibility(:B171, args...; kwargs...)

end # module
