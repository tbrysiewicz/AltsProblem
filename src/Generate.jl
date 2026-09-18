# ============================================================
# 1. Coupler equation G_s(z,w) from the Grassmannian chart
# ============================================================

function generate_G()
    R, t, s, zvec, wvec = polynomial_ring(QQ,
        :t => 0:8, :s => 1:9, :z => 1:1, :w => 1:1)

    z = zvec[1]; w = wvec[1]
    t0,t1,t2,t3,t4,t5,t6,t7,t8 = t
    s1,s2,s3,s4,s5,s6,s7,s8,s9 = s

    M = matrix(R,[
        t0 t3 t6 1 0 0;
        t1 t4 t7 0 1 0;
        t2 t5 t8 0 0 1
    ])

    c00,c01,c02,c11,c12,c22 = [[M[i,j] for i in 1:3] for j in 1:6]
    b00 = [c12[i]*z*w+c01[i]*z+c02[i]*w+c00[i] for i in 1:3]
    b20 = [c22[i]*w+c02[i] for i in 1:3]; b01 = [c11[i]*z+c01[i] for i in 1:3]
    b21 = copy(c12)

    B = matrix(R,hcat(b00,b20,b01,b21))
    minor3(js) = det(B[1:3,js])

    k21 = -minor3([2,3,4]); k01 = minor3([1,3,4])
    k20 = -minor3([1,2,4]); k00 = minor3([1,2,3])

    # Coupler curve F_Lambda(z,w).
    F = k21*k00-k20*k01

    # Pluecker coordinates of the 3-plane represented by M.
    p(i,j,k) = det(M[1:3,[i,j,k]])

    p124 = p(1,2,4); p125 = p(1,2,5); p145 = p(1,4,5)
    p234 = p(2,3,4); p235 = p(2,3,5); p245 = p(2,4,5)
    p246 = p(2,4,6); p256 = p(2,5,6); p345 = p(3,4,5)
    p456 = p(4,5,6)

    graph_relations = [
        p245*s1-p456, p245*s2-p145, p245*s3-p256,
        p245*s4+p246, p245*s5+p235, p245*s6-p234,
        p456*s7-p345, p245*s8+p125, p245*s9-p124
    ]

    Igraph = ideal(R,graph_relations)

    # Restrict to p245*p456 != 0, then eliminate t0,...,t8.
    Iopen = saturation(Igraph+ideal(R,[F]),ideal(R,[p245*p456]))
    Jelim = eliminate(Iopen,collect(t))
    Jgens = [f for f in collect(gens(Jelim)) if !iszero(f)]

    @assert length(Jgens) == 1
    G = -Jgens[1]
    @assert coeff(G,[z,w],[3,3]) == -s1^2

    return (; R,t,s,z,w,G)
end


# ============================================================
# 2. Generic F^(4)(s; lambda)
# ============================================================

function generate_F4()
    dataG = generate_G()
    R = dataG.R; G = dataG.G

    A, ell, ss, u, epsvec, delvec, Zvec, Wvec = polynomial_ring(QQ,
        :ell => (1:9,1:7), :s => 1:9, :u => 1:9,
        :eps => 1:1, :delta => 1:1, :Z => 1:1, :W => 1:1)

    epsilon = epsvec[1]; delta = delvec[1]; Z = Zvec[1]; W = Wvec[1]

    # Map t_i -> 0, s_i -> ss_i, z -> Z, w -> W.
    phi_G = hom(R,A,vcat(fill(zero(A),9),collect(ss),[Z,W]))
    GA = phi_G(G)

    g(i,j) = coeff(GA,[Z,W],[i,j])

    P  = [(0,0),(0,1),(1,0),(0,2),(1,1),(2,0),(0,3),(1,2),(3,0)]
    Pc = [(2,1),(3,1),(2,2),(3,2),(1,3),(2,3),(3,3)]

    Lambda = matrix(A,ell)
    Cp  = matrix(A,9,1,[g(a,b) for (a,b) in P])
    Pcp = matrix(A,7,1,[g(a,b) for (a,b) in Pc])

    F4mat = Cp+Lambda*Pcp
    F4 = [F4mat[i,1] for i in 1:9]

    return (; A,ell,ss,u,epsilon,delta,Z,W,G=GA,P,Pc,F4)
end


# ============================================================
# 3. Generic perturbation G^(4) and homotopy H
# ============================================================

function projected_support(f,variables)
    inds = var_index.(variables)
    unique(Tuple(Int.(a[inds])) for a in exponents(f))
end

monomial_from_exponents(vars,a) =
    prod(vars[i]^a[i] for i in eachindex(vars))

function random_nonzero_integer(rng; bound=9)
    c = 0
    while c == 0
        c = rand(rng,-bound:bound)
    end
    c
end

generic_same_support(f,vars,rng) =
    sum(random_nonzero_integer(rng)*monomial_from_exponents(vars,a)
        for a in projected_support(f,vars))

function generate_homotopy(; seed=20260917)
    data = generate_F4()
    F4 = data.F4; ss = data.ss; epsilon = data.epsilon

    rng = MersenneTwister(seed)
    G4 = [generic_same_support(f,ss,rng) for f in F4]

    expected = [35,38,37,38,43,36,35,41,33]
    @assert [length(projected_support(f,ss)) for f in F4] == expected
    @assert [length(projected_support(f,ss)) for f in G4] == expected

    H4 = [(1-epsilon)*F4[i]+epsilon*G4[i] for i in 1:9]

    return merge(data,(G4=G4,H4=H4,seed=seed))
end