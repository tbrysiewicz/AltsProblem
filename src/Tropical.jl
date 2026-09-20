using HomotopyContinuation

# ============================================================
# Boundary regimes and their compatibility data
# ============================================================

mutable struct BoundaryRegime
    id::Symbol
    name::String
    omega::Vector{Rational{Int}}
    e::Int
    facial_indices::Vector{Int}
    u_indices::Vector{Int}
    substitution::Function
    row_factors::Function       # non-delta factors only
    guard::Function
    count::Int
    source::String
    chart::Union{Nothing,Function}
    K::Any
    orders::Any
    L::Any
    M::Any
    guard_factors::Any
    declared_guard::Function   # stable chart guard; never replaced by a built guard
    S::Any
    coordinate_orders::Any
    construction_metadata::Any
end

# Preserve the existing positional construction API. Runtime construction data
# belong to each regime, while declared_guard always retains the original guard.
function BoundaryRegime(id, name, omega, e, facial_indices, u_indices,
                        substitution, row_factors, guard, count, source, chart,
                        K, orders, L, M, guard_factors)
    BoundaryRegime(id, name, omega, e, facial_indices, u_indices,
        substitution, row_factors, guard, count, source, chart,
        K, orders, L, M, guard_factors, guard, nothing, nothing, nothing)
end

"""
    chart(B::BoundaryRegime)

Return the stored substitution `S_B(u, delta)` over an Oscar fraction field,
with the local variables created automatically.
"""
function chart(B::BoundaryRegime)
    B.chart === nothing && error("$(B.id) has no recorded substitution chart")
    n_u = length(B.omega)
    R, u, delta = polynomial_ring(QQ, :u => 1:n_u, :delta => 1:1)
    F = fraction_field(R)
    B.chart(F.(u), F(delta[1]))
end

function BoundaryRegime(name, omega, e, facial_indices, u_indices,
                        substitution, row_factors, guard)
    BoundaryRegime(Symbol(name), String(name), Rational{Int}.(omega), e,
        facial_indices, u_indices, substitution, row_factors, guard,
        0, "", substitution, nothing, nothing, nothing, nothing, nothing)
end

function BoundaryRegime(id, omega, e, count, source,
                        data::BoundaryRegime)
    label = String(source)
    BoundaryRegime(Symbol(id), data.name, Rational{Int}.(omega), e,
        data.facial_indices, data.u_indices, data.substitution,
        data.row_factors, data.declared_guard, count, label, data.substitution,
        nothing, nothing, nothing, nothing, nothing)
end

function BoundaryRegime(id, omega, e, count, source,
                        chart::Function)
    label = String(source)
    BoundaryRegime(Symbol(id), String(id), Rational{Int}.(omega), e,
        Int[], Int[], chart, _no_extra_factors, X -> one(X[1]),
        count, label, chart, nothing, nothing, nothing, nothing, nothing)
end

function BoundaryRegime(id, omega, e, count, source,
                        ::Nothing, chart::Function)
    BoundaryRegime(id, omega, e, count, source, chart)
end

function BoundaryRegime(id, omega, e, count, source)
    label = String(source)
    BoundaryRegime(Symbol(id), String(id), Rational{Int}.(omega), e,
        Int[], Int[], X -> error("no chart recorded for $(id)"),
        _no_extra_factors, X -> one(X[1]), count, label, nothing,
        nothing, nothing, nothing, nothing, nothing)
end

function Base.show(io::IO,::MIME"text/plain",B::BoundaryRegime)
    println(io,"Boundary regime: ",B.name)
    println(io,"  omega          = ",B.omega)
    println(io,"  e              = ",B.e)
    println(io,"  facial         = u",B.facial_indices)
    println(io,"  transverse/etc = u",B.u_indices)

    D,uu,dd = polynomial_ring(QQ,:u => 1:9,:delta => 1:1)
    K = fraction_field(D); uK = K.(uu); deltaF = K(dd[1])

    if B.guard_factors === nothing
        println(io,"  guard          = ",B.guard(uK))
    else
        println(io,"  guard factors  = ",B.guard_factors)
    end

    println(io,"\n  substitution:")
    for (i,f) in enumerate(B.substitution(uK,deltaF))
        println(io,"    s[$i] = ",f)
    end

    factors = B.row_factors(uK,deltaF)
    println(io,"\n  extra row factors:")
    if all(isone,factors)
        println(io,"    none")
    elseif all(f -> f == factors[1],factors)
        println(io,"    all rows: ",factors[1])
    else
        for (i,f) in enumerate(factors)
            println(io,"    H[$i]: ",f)
        end
    end
end

Base.getproperty(B::BoundaryRegime, name::Symbol) =
    name === :ramification ? getfield(B, :e) : getfield(B, name)


# ============================================================
# Finite compatibility systems
# ============================================================

function delta_order(f,delta)
    iszero(f) && error("An identically zero row needs a different compatibility construction")
    i = var_index(delta)
    ord(p) = minimum(Int(a[i]) for a in exponents(p))
    ord(numerator(f))-ord(denominator(f))
end

function construct_scaled_boundary(B::BoundaryRegime, H, s, u, epsilon, delta, scales;
                                   automatic_orders=false)
    A = parent(H[1]); Q = fraction_field(A)
    uq, dq = Q.(u), Q(delta)
    row_scales = scales(uq, dq)
    length(row_scales) == length(H) || error("$(B.id) needs one scale per row")
    images = Q.(gens(A))
    images[var_index.(s)] = B.chart(uq, dq)
    images[var_index(epsilon)] = dq^B.e
    phi = hom(A, Q, images)
    substituted = phi.(H)
    initial_orders = [delta_order(f, delta) for f in substituted]
    rows = [row_scales[i]*substituted[i] for i in eachindex(H)]
    orders = automatic_orders ? [delta_order(f, delta) for f in rows] :
                               zeros(Int, length(H))
    M = zero_matrix(Q, length(H), length(H))
    L = similar(rows)
    for i in eachindex(H)
        M[i,i] = row_scales[i]/dq^orders[i]
        L[i] = rows[i]/dq^orders[i]
    end
    out = (; L, M, row_orders=orders, initial_orders,
             guard_factors=typeof(zero(A))[])
    return finish_boundary_construction!(B, out, H, s, u, epsilon, delta)
end

function finite_compatibility_system(B::BoundaryRegime, H, s, u, epsilon, delta)
    _, _, K, _ = normalize_boundary(B, H, s, u, epsilon, delta)
    return K, ideal(parent(H[1]), K), B.orders
end


# ============================================================
# Oscar -> HomotopyContinuation.jl
# ============================================================

function toHC(f,oscarvars,hcvars)
    length(oscarvars) == length(hcvars) ||
        throw(ArgumentError("Oscar/HC symbol counts differ"))
    inds = var_index.(oscarvars)
    length(unique(inds)) == length(inds) ||
        throw(ArgumentError("Each Oscar symbol must be bound exactly once"))
    bound = Set(inds)
    for a in exponents(f)
        any(a[j] != 0 for j in eachindex(a) if !(j in bound)) &&
            throw(ArgumentError("Unbound symbol in Oscar-to-HC conversion"))
    end
    iszero(f) && return 0
    sum(begin
        c = coeff(f,i); a = exponent_vector(f,i)
        cq = Rational{BigInt}(BigInt(numerator(c)),BigInt(denominator(c)))
        cq*prod(hcvars[j]^Int(a[inds[j]]) for j in eachindex(hcvars))
    end for i in 1:length(f))
end
function active_variables(K)
    R = parent(first(K)); gs = collect(gens(R))
    used = falses(length(gs))
    for f in K, a in exponents(f)
        used .|= a .!= 0
    end
    gs[used]
end

function solve_compatibility_system(K,parameters; seed=20260917)
    rng = MersenneTwister(seed)
    vals = [random_nonzero_integer(rng) for _ in parameters]
    Kspec = [Oscar.evaluate(f,parameters,vals) for f in K]

    variables = active_variables(Kspec)
    hcvars = [Variable(:u,i) for i in eachindex(variables)]

    Khc = [toHC(f,variables,hcvars) for f in Kspec]
    Fhc = System(Khc; variables=hcvars)

    return HomotopyContinuation.solve(Fhc),vals,variables,hcvars
end

function hc_compatibility_system(B::BoundaryRegime; keep_parameters=true, seed=20260917)
    data = generate_homotopy(; seed=seed)

    K, J, orders = finite_compatibility_system(
        B,
        data.H4,
        data.ss,
        data.u,
        data.epsilon,
        data.delta,
    )

    return hc_compatibility_system(K,vec(data.ell);
        keep_parameters=keep_parameters, seed=seed, variables=data.u)
end

function hc_compatibility_system(K, parameters; keep_parameters=false,
                                 seed=20260917, variables=nothing)
    parameters = collect(parameters)
    unknowns = isnothing(variables) ?
        [v for v in active_variables(K) if !(v in parameters)] : collect(variables)
    length(unknowns) == length(K) ||
        throw(ArgumentError("Compatibility system must have one equation per local variable"))
    isempty(intersect(unknowns, parameters)) ||
        throw(ArgumentError("Local variables and parameters must be disjoint"))
    hcvars = HomotopyContinuation.ModelKit.variables(:u, 1:length(unknowns))
    if keep_parameters
        hcparams = HomotopyContinuation.ModelKit.variables(:l, 1:length(parameters))
        Khc = [toHC(f, vcat(unknowns, parameters), vcat(hcvars, hcparams)) for f in K]
        return System(Khc; variables=hcvars, parameters=hcparams)
    end
    rng = MersenneTwister(seed)
    vals = [random_nonzero_integer(rng) for _ in parameters]
    Kspec = [Oscar.evaluate(f, parameters, vals) for f in K]
    Khc = [toHC(f, unknowns, hcvars) for f in Kspec]
    return System(Khc; variables=hcvars)
end
