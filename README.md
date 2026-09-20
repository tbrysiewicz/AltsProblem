# AltCertify

`AltCertify` certifies contributions from boundary regimes to the BKK defect of
the formulation of the nine-point path-synthesis problem in Alt's problem. It
also certifies a lower bound of 8,652 mechanisms, and hence 1,442 coupler-curve
nets, through nine prescribed points.

## Setup

Start Julia from the repository root and activate the project:

```julia-repl
$ julia

pkg> activate .
pkg> instantiate

julia> using AltCertify
```

The first run may take longer while Julia downloads dependencies and compiles
the package.

## Inspect a boundary regime

Extract one of the boundary regimes:

```julia-repl
julia> B = B11
Boundary regime: B_1, regime 1
omega
= Rational{Int64}[1//2, 0, 0, 0, 0, 0, -1//2, 0, 1//2]
e
= 2
facial
= u[1, 2, 5, 6]
transverse/etc = u[3, 4, 7, 8, 9]
guard
= u[1]*u[2]*u[5]*u[6]
substitution:
s[1] = u[1]*delta[1]
s[2] = u[4]*delta[1] - 2*u[6]
s[3] = -u[1]*u[6]*delta[1] - u[3]*u[5]*delta[1] + u[5]^2 + u[8]*delta[1]^2
s[4] = u[3]*delta[1] - 2*u[5]
s[5] = -u[5]*u[6] + u[7]*delta[1]^2
s[6] = u[6]
s[7] = u[5]//(u[1]*delta[1])
s[8] = -u[2]*u[5]*delta[1] + u[4]*u[6]*delta[1] - u[6]^2 + u[9]*delta[1]^2
s[9] = u[2]*delta[1]
extra row factors:
none
```

## Certify a boundary contribution

Certify the contribution of this boundary regime to the BKK defect:

```julia-repl
julia> certify(B11)

1. MONODROMY SOLVE
   monodromy solutions: 1322
2. SOFT CERTIFICATION AT FLOAT64 PARAMETERS
   soft-certified distinct solutions: 1322
3. NEARBY EXACT RATIONAL PARAMETERS
   maximum parameter displacement: 1.4060347687259386e-12
   exact parameter vector constructed.
4. PARAMETER HOMOTOPY TO RATIONAL FIBER
   finite endpoint candidates for hard certification: 1322
5. HARD CERTIFICATION AT EXACT RATIONAL PARAMETERS
   hard-certified distinct solutions: 1322
6. GUARD CHECK ON CERTIFIED BOXES
   guard excludes zero on 1322 / 1322 boxes
   PASS: 1322 distinct exact-rational certified solutions, all guarded.
```

By default, the exact system, construction metadata, and interval certificates
are written beneath `Certificates/B11`.

Every boundary follows this same API. `normalize_boundary(B)` returns
`(L, M, K, guard)` and stores its substitution as `B.S`. The identities are
`L = M*H(S,delta^B.e)` and `K = L(u,0)`. After construction, `guard` returns a
list of factor values; parameter-dependent guards also accept the original
parameter symbols and their exact values. `certify(B)` binds those parameters
automatically. See [the boundary audit](BOUNDARY_AUDIT.md) for the necessary
differences in row operations and the verification commands.

`reproduce_certification(B)` performs a new solve using the saved seed and
overwrites that boundary's output after success; it does not replay the old
interval certificates.

## Certify the lower bound for Alt's problem

Run the full certification to prove the existence of at least 8,652 mechanisms
drawing coupler curves through nine points, and therefore at least 1,442 nets:

```julia-repl
julia> certify_alt_mechanisms()
Finding 8652 lower-bound solutions by monodromy...
Solutions found: 8652    Time: 0:00:40
tracked loops (queued): 17156 (132)
solutions in current (last) loop: 8638 (13)
generated loops (no change): 2 (0)
Tracking solutions to exact rational parameters...
Tracking 8652 paths... 100%|████████████████████████████████████████| Time: 0:00:08
# paths tracked: 8652

# non-singular solutions (real): 8652 (2124)
# singular endpoints (real): 0 (0)
# total solutions (real): 8652 (2124)

Certifying 8652 tracked solutions...
Certifying 8652 solutions... 100%|████████████████████████████████████████| Time: 0:00:01
# processed: 8652

# certified (real): 8652 (2124)
# distinct (real): 8652 (2124)

Writing certificates
Saved 8652 distinct certificates.
CertificationResult
===================

• 8652 solution candidates given
• 8652 certified solution intervals (2124 real, 6528 complex)
• 8652 distinct certified solution intervals (2124 real, 6528 complex)
```

Progress, timings, and the displayed certificate path vary by machine.
