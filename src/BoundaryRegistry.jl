const BOUNDARY_REGIMES = [
    BoundaryRegime(:B11, boundary_B11.omega, 2, 1322, "BoundaryData/B11", boundary_B11),
    BoundaryRegime(:B21, boundary_B21.omega, 1, 335, "BoundaryData/B21", boundary_B21),
    BoundaryRegime(:B31, boundary_B31.omega, 1, 1648, "BoundaryData/B31", boundary_B31),
    BoundaryRegime(:B32, boundary_B32.omega, 2, 542, "BoundaryData/B32", boundary_B32),
    BoundaryRegime(:B41, boundary_B41.omega, 1, 6, "BoundaryData/B41", boundary_B41),
    BoundaryRegime(:B42, boundary_B42.omega, 2, 14, "BoundaryData/B42", boundary_B42),
    BoundaryRegime(:B43, boundary_B43.omega, 3, 54, "BoundaryData/B43", boundary_B43),
    BoundaryRegime(:B44, boundary_B44.omega, 5, 30, "BoundaryData/B44", boundary_B44),
    BoundaryRegime(:B51, boundary_B51.omega, 6, 18, "BoundaryData/B51", boundary_B51),
    BoundaryRegime(:B52, boundary_B52.omega, 10, 20, "BoundaryData/B52", boundary_B52),
    BoundaryRegime(:B53, boundary_B53.omega, 6, 12, "BoundaryData/B53", boundary_B53),
    BoundaryRegime(:B61, boundary_B61.omega, 4, 16, "BoundaryData/B61", boundary_B61),
    BoundaryRegime(:B71, boundary_B71.omega, 4, 16, "BoundaryData/B71", boundary_B71),
    BoundaryRegime(:B81, boundary_B81.omega, 2, 4, "BoundaryData/B81", boundary_B81),
    BoundaryRegime(:B91, [1//1,-1//1,2//1,1//1,2//1,1//1,0//1,0//1,-1//1], 1, 8, "BoundaryData/B91", nothing, S_B91),
    BoundaryRegime(:B101, boundary_B101.omega, 4, 8, "BoundaryData/B101", boundary_B101),
    BoundaryRegime(:B111, boundary_B111.omega, 6, 6, "BoundaryData/B111", boundary_B111),
    BoundaryRegime(:B112, boundary_B112.omega, 10, 10, "BoundaryData/B112", boundary_B112),
    BoundaryRegime(:B121, boundary_B121.omega, 3, 6, "BoundaryData/B121", boundary_B121),
    BoundaryRegime(:B122, boundary_B122.omega, 4, 4, "BoundaryData/B122", boundary_B122),
    BoundaryRegime(:B131, [-1//1,-2//3,-2//1,-1//1,-5//3,-2//3,0//1,-4//3,-1//3], 3, 6, "BoundaryData/B131", nothing, S_B131),
    BoundaryRegime(:B141, [-1//1,-3//4,-2//1,-1//1,-7//4,-3//4,0//1,-3//2,-1//2], 4, 4, "BoundaryData/B141", nothing, S_B141),
    BoundaryRegime(:B151, [-4//1,-5//2,-5//1,-3//1,-3//1,-1//1,2//1,-7//2,-3//2], 2, 2, "BoundaryData/Boundary15To18", nothing, S_B151),
    BoundaryRegime(:B161, [-2//1,-1//1,-3//1,-2//1,-2//1,-1//1,1//1,-2//1,-1//1], 1, 2, "BoundaryData/Boundary15To18", nothing, S_B161),
    BoundaryRegime(:B171, [-4//1,-2//1,-5//1,-4//1,-2//1,-1//1,3//1,-3//1,-2//1], 1, 1, "BoundaryData/Boundary15To18", nothing, S_B171),
    BoundaryRegime(:B181, [-3//2,-1//2,-3//1,-2//1,-5//2,-3//2,1//2,-2//1,-1//1], 2, 2, "BoundaryData/Boundary15To18", nothing, S_B181),
]

boundary_regimes() = BOUNDARY_REGIMES

function boundary_regime(id::Symbol)
    i = findfirst(b -> b.id == id, BOUNDARY_REGIMES)
    isnothing(i) && throw(KeyError(id))
    BOUNDARY_REGIMES[i]
end

@assert length(BOUNDARY_REGIMES) == 26
@assert sum(b.count for b in BOUNDARY_REGIMES) == 4096
@assert length(unique(b.id for b in BOUNDARY_REGIMES)) == 26

# Public REPL names. Every B... binding has the same BoundaryRegime type.
const B11 = boundary_regime(:B11)
const B21 = boundary_regime(:B21)
const B31 = boundary_regime(:B31)
const B32 = boundary_regime(:B32)
const B41 = boundary_regime(:B41)
const B42 = boundary_regime(:B42)
const B43 = boundary_regime(:B43)
const B44 = boundary_regime(:B44)
const B51 = boundary_regime(:B51)
const B52 = boundary_regime(:B52)
const B53 = boundary_regime(:B53)
const B61 = boundary_regime(:B61)
const B71 = boundary_regime(:B71)
const B81 = boundary_regime(:B81)
const B91 = boundary_regime(:B91)
const B101 = boundary_regime(:B101)
const B111 = boundary_regime(:B111)
const B112 = boundary_regime(:B112)
const B121 = boundary_regime(:B121)
const B122 = boundary_regime(:B122)
const B131 = boundary_regime(:B131)
const B141 = boundary_regime(:B141)
const B151 = boundary_regime(:B151)
const B161 = boundary_regime(:B161)
const B171 = boundary_regime(:B171)
const B181 = boundary_regime(:B181)
