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
    label = "BoundaryData/" * String(id)
    BoundaryRegime(Symbol(id), data.name, Rational{Int}.(omega), e,
        data.facial_indices, data.u_indices, data.substitution,
        data.row_factors, data.guard, count, label, data.substitution,
        nothing, nothing, nothing, nothing, nothing)
end

function BoundaryRegime(id, omega, e, count, source,
                        chart::Function)
    label = "BoundaryData/" * String(id)
    BoundaryRegime(Symbol(id), String(id), Rational{Int}.(omega), e,
        Int[], Int[], chart, _no_extra_factors, X -> one(X[1]),
        count, label, chart, nothing, nothing, nothing, nothing, nothing)
end

function BoundaryRegime(id, omega, e, count, source,
                        ::Nothing, chart::Function)
    BoundaryRegime(id, omega, e, count, source, chart)
end

function BoundaryRegime(id, omega, e, count, source)
    label = "BoundaryData/" * String(id)
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
    i = var_index(delta)
    ord(p) = minimum(Int(a[i]) for a in exponents(p))
    ord(numerator(f))-ord(denominator(f))
end

function construct_scaled_boundary(B::BoundaryRegime, H, s, u, epsilon, delta, scales)
    A = parent(H[1])
    Q = fraction_field(A)
    uQ = Q.(u)
    deltaQ = Q(delta)
    row_scales = scales(uQ, deltaQ)
    length(row_scales) == length(H) || error("$(B.id) needs one scale per row")

    images = Q.(gens(A))
    images[var_index.(s)] = B.chart(uQ, deltaQ)
    images[var_index(epsilon)] = deltaQ^B.e
    phi = hom(A, Q, images)

    M = zero_matrix(Q, length(H), length(H))
    L = Vector{typeof(deltaQ)}(undef, length(H))
    K = Vector{typeof(zero(A))}(undef, length(H))
    for i in eachindex(H)
        M[i, i] = row_scales[i]
        row = row_scales[i] * phi(H[i])
        common = gcd(numerator(row), denominator(row))
        row = Q(divexact(numerator(row), common)) /
              Q(divexact(denominator(row), common))
        den0 = Oscar.evaluate(denominator(row), [delta], [zero(A)])
        iszero(den0) && error("$(B.id) row $i is not regular at delta=0")
        num0 = Oscar.evaluate(numerator(row), [delta], [zero(A)])
        K[i] = numerator(Q(num0) / Q(den0))
        L[i] = row
    end

    B.L = L
    B.M = M
    B.K = K
    B.orders = fill(0, length(H))
    return L, M, K, B.guard
end

function finite_compatibility_system(B::BoundaryRegime,H,s,u,epsilon,delta)
    if B.id === :B42
        _, _, K, _ = construct_B42(B, H, s, u, epsilon, delta)
        return K, ideal(parent(H[1]), K), B.orders
    elseif B.id === :B43
        _, _, K, _ = construct_B43(B, H, s, u, epsilon, delta)
        return K, ideal(parent(H[1]), K), B.orders
    elseif B.id === :B44
        _, _, K, _ = construct_B44(B, H, s, u, epsilon, delta)
        return K, ideal(parent(H[1]), K), B.orders
    elseif B.id in RECORDED_COMPATIBILITY_IDS
        _, _, K, _ = construct_recorded_boundary(B, H, s, u, epsilon, delta)
        return K, ideal(parent(H[1]), K), B.orders
    elseif B.id in GENERAL_COMPATIBILITY_IDS
        _, _, K, _ = construct_general_boundary(B, H, s, u, epsilon, delta)
        return K, ideal(parent(H[1]), K), B.orders
    end
    isempty(B.facial_indices) &&
        error("$(B.id) has no compatibility constructor")
    A = parent(H[1]); Kfrac = fraction_field(A)
    uK = Kfrac.(u); deltaF = Kfrac(delta)

    images = Kfrac.(gens(A))
    images[var_index.(s)] = B.substitution(uK,deltaF)
    images[var_index(epsilon)] = deltaF^B.e
    phi = hom(A,Kfrac,images)

    Hsub = phi.(H)

    # Remove boundary-specific non-delta factors first.
    extras = B.row_factors(uK,deltaF)
    Hsub = [Hsub[i]/extras[i] for i in eachindex(Hsub)]

    # Find and remove the lowest delta power automatically.
    orders = [delta_order(h,delta) for h in Hsub]
    Hscaled = [Hsub[i]/deltaF^orders[i] for i in eachindex(Hsub)]

    # Set delta = 0.
    K = map(Hscaled) do f
        num = Oscar.evaluate(numerator(f),[delta],[zero(A)])
        den = Oscar.evaluate(denominator(f),[delta],[zero(A)])
        @assert !iszero(den)
        numerator(Kfrac(num)/Kfrac(den))
    end

    B.K = K
    B.orders = orders
    return K, ideal(A,K), orders
end


# ============================================================
# Oscar -> HomotopyContinuation.jl
# ============================================================

function toHC(f,oscarvars,hcvars)
    inds = var_index.(oscarvars)
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
        keep_parameters=keep_parameters, seed=seed)
end

function hc_compatibility_system(K,parameters; keep_parameters=false, seed=20260917)
    R = parent(first(K))
    gensR = collect(gens(R))
    parameter_positions = Set(var_index.(parameters))

    if keep_parameters
        # Find the active nonparameter Oscar variables.
        used = falses(length(gensR))
        for f in K, a in exponents(f)
            used .|= a .!= 0
        end

        variable_positions = [i for i in eachindex(gensR)
            if used[i] && !(i in parameter_positions)]

        variables = gensR[variable_positions]

        # Make separate HC variables and parameters.
        hcvars = HomotopyContinuation.ModelKit.variables(:u,1:length(variables))
        hcparams = HomotopyContinuation.ModelKit.variables(:l,1:length(parameters))

        oscar_symbols = vcat(variables,collect(parameters))
        hc_symbols = vcat(hcvars,hcparams)

        Khc = [toHC(f,oscar_symbols,hc_symbols) for f in K]

        return System(Khc; variables=hcvars, parameters=hcparams)

    else
        rng = MersenneTwister(seed)
        vals = [random_nonzero_integer(rng) for _ in parameters]
        Kspec = [Oscar.evaluate(f,parameters,vals) for f in K]

        used = falses(length(gensR))
        for f in Kspec, a in exponents(f)
            used .|= a .!= 0
        end

        variables = gensR[used]
        hcvars = HomotopyContinuation.ModelKit.variables(:u,1:length(variables))

        Khc = [toHC(f,variables,hcvars) for f in Kspec]

        return System(Khc; variables=hcvars)
    end
end
