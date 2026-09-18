module AltCertify

using Oscar
using Random
using HomotopyContinuation

include("Generate.jl")
include("Tropical.jl")
include("Boundaries.jl")

include("BoundaryData/B11.jl")
include("BoundaryData/B21.jl")
include("BoundaryData/B31.jl")
include("BoundaryData/B32.jl")
include("BoundaryData/B41.jl")
include("BoundaryData/B42.jl")
include("BoundaryData/B43.jl")
include("BoundaryData/B44.jl")
include("BoundaryData/B51.jl")
include("BoundaryData/B52.jl")
include("BoundaryData/B53.jl")
include("BoundaryData/B61.jl")
include("BoundaryData/B71.jl")
include("BoundaryData/B81.jl")
include("BoundaryData/B91.jl")
include("BoundaryData/B101.jl")
include("BoundaryData/B111.jl")
include("BoundaryData/B112.jl")
include("BoundaryData/B121.jl")
include("BoundaryData/B122.jl")
include("BoundaryData/B131.jl")
include("BoundaryData/B141.jl")
include("BoundaryData/Boundary15To18.jl")

include("BoundaryCompatibilityExtras.jl")
include("BoundaryCompatibilityGeneral.jl")
include("BoundaryRegistry.jl")
include("BoundaryCompatibility.jl")
include("normalize_boundary.jl")
include("Certify.jl")

using .Certify: certify_regime, CertificationRun, guard_excludes_zero

export generate_G,generate_F4,generate_homotopy
export finite_compatibility_system
export normalize_boundary
export solve_compatibility_system
export BoundaryRegime, BOUNDARY_REGIMES, boundary_regime, boundary_regimes
export chart
export B11, B21, B31, B32, B41, B42, B43, B44
export B51, B52, B53, B61, B71, B81, B91, B101
export B111, B112, B121, B122, B131, B141, B151, B161, B171, B181
export hc_compatibility_system
export certify_regime, CertificationRun, guard_excludes_zero
using .Certify: certify, reproduce_certification
export certify, reproduce_certification

end
