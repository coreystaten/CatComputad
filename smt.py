from fill import *
from ontology import *

from unitor import *
from shuffle import *
from syllepsis import *
from tricathom import *

_H = Functor("H")

_a = ConstPrim0("a")
_b = ConstPrim0("b")
_c = ConstPrim0("c")
_d = ConstPrim0("d")
_e = ConstPrim0("e")

_Ha = app(_H, _a)
_Hb = app(_H, _b)
_Hc = app(_H, _c)
_Hd = app(_H, _d)
_He = app(_H, _e)


# TODO: Choose to use _H1x instead of _Hx whenever possible, since it's easier to reduce than expand?
_H1a = app(_H, ensureEqMol1(_a))
_H1b = app(_H, ensureEqMol1(_b))
_H1c = app(_H, ensureEqMol1(_c))
_H1d = app(_H, ensureEqMol1(_d))
_H1e = app(_H, ensureEqMol1(_e))

chi1Source = minimalASTFromMol0(ensureMol0(
    comp0(app(_H,_a), app(_H, _b))), [_a, _b], [_H])
chi1Target = minimalASTFromMol0(ensureMol0(
    app(_H, comp0(_a, _b))), [_a, _b], [_H])
chi1 = PrimitiveFamily("chi1", 1, 1, [0, 0], chi1Source, chi1Target)
chi1Dim2 = tritrans2CellPrimFamily(chi1)
chi1Dim3 = tritrans3CellPrimFamily(chi1, chi1Dim2)
chi1Pi = tritransPiPrimFamily(chi1, chi1Dim2)
chi1M = tritransMPrimFamily(chi1, chi1Dim2)

# We use the same prim for the unit in both categories; this shouldn't cause any problems, since these we can tell which is meant by which category it should fall in.
iota1Source = ConstNode(unit)
iota1Target = FunctorNode(0, ConstNode(unit))
iota1 = PrimitiveFamily("iota1", 1, 1, [], iota1Source, iota1Target)

omega2Source = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(_Ha, chi1.fprimf([_H], _b, _c)),
        chi1.fprimf([_H], _a, comp0(_b, _c)))), [_a, _b, _c], [_H])
omega2Target = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(chi1.fprimf([_H], _a, _b), _Hc),
        chi1.fprimf([_H], comp0(_a, _b), _c))), [_a, _b, _c], [_H])
omega2 = PrimitiveFamily("omega2", 2, 1, [0, 0, 0], omega2Source, omega2Target)
omega2Dim3 = trimod3CellPrimFamily(omega2, {chi1: chi1Dim2})


gamma2Source = minimalASTFromEqMol1(ensureEqMol1(
    comp1s(
        comp0(iota1.fprimf([_H]), app(_H, _a)),
        chi1.fprimf([_H], unit, _a),
        app(_H, unitor10.fprim(_a))
    )), [_H], [_a])
gamma2Target = minimalASTFromEqMol1(ensureEqMol1(
    unitor10.fprim(app(_H, _a))), [_a], [_H])
gamma2 = PrimitiveFamily("gamma2", 2, 1, [0], gamma2Source, gamma2Target)

delta2Source = minimalASTFromEqMol1(ensureEqMol1(
    comp1s(
        unitor11.adj.fprim(app(_H, _a)),
        comp0(app(_H, _a), iota1.fprimf([_H])),
        chi1.fprimf([_H], _a, unit)
    )), [_a], [_H])
delta2Target = minimalASTFromEqMol1(ensureEqMol1(
    app(_H, unitor11.adj.fprim(_a))), [_a], [_H])
delta2 = PrimitiveFamily("delta2", 2, 1, [0], delta2Source, delta2Target)

# We switch a and b form other notes.
# We also switch source and target, since it's u2.adj that seems to come up in diagrams otherwise.
u2Source = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        shuffle1.fprim(app(_H, _a), app(_H, _b)),
        chi1.fprimf([_H], _b, _a))), [_a, _b], [_H])
u2Target = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        chi1.fprimf([_H], _a, _b),
        app(_H, shuffle1.fprim(_a, _b)))), [_a, _b], [_H])
u2 = PrimitiveFamily("u2", 2, 1, [0, 0], u2Source, u2Target)
u2Dim3 = trimod3CellPrimFamily(u2, {chi1: chi1Dim2, shuffle1: shuffle1Dim2})

#omega3Source0 = ensureMol0(comp0s(app(_H, _a), app(_H, _b), app(_H, _c), app(_H, _d)))
#omega3Target0 = ensureMol0(app(_H, comp0s(_a, _b, _c, _d)))
#omega3Paths1 = findPaths1(omega3Source0, omega3Target0, [chi1])
#searchForPathPairs2(omega3Paths1, [omega2])
#[((H(a) @ H(b)) @ chi1{0}(c, d)) . (H(a) @ chi1{0}(b, (c @ d))) . (chi1{0}(a, (b @ c @ d)))]
# --
#[(chi1{H}(a, b) @ (H(c) @ H(d))) . (chi1{H}((a @ b), c) @ H(d)) . (chi1{0}((a @ b @ c), d))]
omega3Source1 = ensureEqMol1(
    comp1s(
        comp0s(_H1a, _H1b, chi1.fprimf([_H], _c, _d)),
        comp0(_H1a, chi1.fprimf([_H], _b, comp0(_c, _d))),
        chi1.fprimf([_H], _a, comp0s(_b, _c, _d))))
omega3Target1 = ensureEqMol1(
    comp1s(
        comp0s(chi1.fprimf([_H], _a, _b), _H1c, _H1d),
        comp0(chi1.fprimf([_H], comp0(_a, _b), _c), _H1d),
        chi1.fprimf([_H], comp0s(_a, _b, _c), _d)))
#omega3Paths2 = findPaths2(omega3Source1, omega3Target1, [omega2])
print(len(cellsAway2(omega3Source1, [omega2])))

# [[((H(a) @ H(b)) @ chi1{H}(c, d)) . omega2{0}(a, b, (c @ d))] &
# [(chi1{H}(a, b) @ (H(c) @ H(d))) . omega2{0}((a @ b), c, d)]]
# --
# [[(H(a) @ omega2{0}(b, c, d)) . (chi1{H}(a, (b @ c @ d)))] &
# [(H(a) @ chi1{H}(b, c) @ H(d)) . omega2{0}(a, (b @ c), d)] &
# [(omega2{0}(a, b, c) @ H(d)) . (chi1{H}((a @ b @ c), d))]]
omega3SourceAST = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        comp1(
            comp0s(app(_H, _a), app(_H, _b), chi1.fprimf([_H], _c, _d)),
            omega2.fprimf([_H], _a, _b, comp0(_c, _d))),
        comp1(
            comp0s(chi1.fprimf([_H], _a, _b), app(_H, _c), app(_H, _d)),
            omega2.fprimf([_H], comp0(_a, _b), _c, _d))
    )
), [_a, _b, _c, _d], [_H])
omega3TargetAST = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1(
            comp0(app(_H, _a), omega2.fprimf([_H], _b, _c, _d)),
            chi1.fprimf([_H], _a, comp0s(_b, _c, _d))),
        comp1(
            comp0s(app(_H, _a), chi1.fprimf([_H], _b, _c), app(_H, _d)),
            omega2.fprimf([_H], _a, comp0(_b, _c), _d)),
        comp1(
            comp0(omega2.fprimf([_H], _a, _b, _c), app(_H, _d)),
            chi1.fprimf([_H], comp0s(_a, _b, _c), _d))
    )
), [_a, _b, _c, _d], [_H])
omega3 = PrimitiveFamily("omega3", 3, 1, [0, 0, 0, 0], omega3SourceAST, omega3TargetAST)

omega4Source0 = ensureMol0(comp0s(app(_H, _a), app(_H, _b), app(_H, _c), app(_H, _d), app(_H, _e)))
omega4Target0 = ensureMol0(app(_H, comp0s(_a, _b, _c, _d, _e)))
#omega4Paths1 = findPaths1(omega4Source0, omega4Target0, [chi1])
#searchForPathPairs2(omega4Paths1, [omega2])
# We pick the choice that lines up with the pattern one dimension down; grouping left versus grouping right.
# This is also the choice with the most 2-dimensional paths; in general, the ability to find properly oriented 3-paths depends on choice of both 1st and 2nd dimension?
# TODO: Do we need to revisit 1-paths for permutahedron?
# [((H(a) @ H(b) @ H(c)) @ chi1{0}(d, e)) . ((H(a) @ H(b)) @ chi1{0}(c, (d @ e))) . (H(a) @ chi1{0}(b, (c @ d @ e))) . (chi1{0}(a, (b @ c @ d @ e)))]
# --
#[(chi1{H}(a, b) @ (H(c) @ H(d) @ H(e))) . (chi1{H}((a @ b), c) @ (H(d) @ H(e))) . (chi1{H}((a @ b @ c), d) @ H(e)) . (chi1{0}((a @ b @ c @ d), e))]
omega4Source1 = ensureEqMol1(
    comp1s(
        comp0s(app(_H, _a), app(_H, _b), app(_H, _c), chi1.fprimf([_H], _d, _e)),
        comp0s(app(_H, _a), app(_H, _b), chi1.fprimf([_H], _c, comp0(_d, _e))),
        comp0s(app(_H, _a), chi1.fprimf([_H], _b, comp0s(_c, _d, _e))),
        chi1.fprimf([_H], _a, comp0s(_b, _c, _d, _e))
    )
)
omega4Target1 = ensureEqMol1(
    comp1s(
        comp0s(chi1.fprimf([_H], _a, _b), app(_H, _c), app(_H, _d), app(_H, _e)),
        comp0s(chi1.fprimf([_H], comp0(_a, _b), _c), app(_H, _d), app(_H, _e)),
        comp0s(chi1.fprimf([_H], comp0s(_a, _b, _c), _d), app(_H, _e)),
        chi1.fprimf([_H], comp0s(_a, _b, _c, _d), _e)
    )
)
#omega4Paths2 = findPaths2(omega4Source1, omega4Target1, [omega2])
#searchForPathPairs3(omega4Paths2, omega4Paths2, [omega3])
# SUCCESS!!!!!!
# [[((H(a) @ H(b) @ H(c)) @ chi1{H}(d, e)) . ((H(a) @ H(b)) @ chi1{H}(c, (d @ e))) . omega2{0}(a, b, (c @ d @ e))] &
# [((H(a) @ H(b) @ H(c)) @ chi1{H}(d, e)) . (chi1{H}(a, b) @ (H(c) @ H((d @ e)))) . omega2{0}((a @ b), c, (d @ e))] &
# [(chi1{H}(a, b) @ (H(c) @ H(d) @ H(e))) . (chi1{H}((a @ b), c) @ (H(d) @ H(e))) . omega2{0}((a @ b @ c), d, e)]]
omega4Source2 = ensureEqAEMol2(
    comp2s(
        comp1s(
            comp0s(app(_H, _a), app(_H, _b), app(_H, _c), chi1.fprimf([_H], _d, _e)),
            comp0s(app(_H, _a), app(_H, _b), chi1.fprimf([_H], _c, comp0(_d, _e))),
            omega2.fprimf([_H], _a, _b, comp0s(_c, _d, _e))),
        comp1s(
            comp0s(app(_H, _a), app(_H, _b), app(_H, _c), chi1.fprimf([_H], _d, _e)),
            comp0s(chi1.fprimf([_H], _a, _b), app(_H, _c), app(_H, comp0(_d, _e))),
            omega2.fprimf([_H], comp0(_a, _b), _c, comp0(_d, _e))),
        comp1s(
            comp0s(chi1.fprimf([_H], _a, _b), app(_H, _c), app(_H, _d), app(_H, _e)),
            comp0s(chi1.fprimf([_H], comp0(_a, _b), _c), app(_H, _d), app(_H, _e)),
            omega2.fprimf([_H], comp0s(_a, _b, _c), _d, _e))
    )
)

# --
# [[((H(a) @ H(b)) @ omega2{0}(c, d, e)) . (H(a) @ chi1{H}(b, (c @ d @ e))) . (chi1{H}(a, (b @ c @ d @ e)))] &
# [((H(a) @ H(b)) @ chi1{H}(c, d) @ H(e)) . (H(a) @ omega2{0}(b, (c @ d), e)) . (chi1{H}(a, (b @ c @ d @ e)))] &
# [(H(a) @ omega2{0}(b, c, d) @ H(e)) . (H(a) @ chi1{0}((b @ c @ d), e)) . (chi1{H}(a, (b @ c @ d @ e)))] &
# [(H(a) @ chi1{H}(b, c) @ (H(d) @ H(e))) . (H(a) @ chi1{0}((b @ c), d) @ H(e)) . omega2{0}(a, (b @ c @ d), e)] &
# [(H(a) @ chi1{H}(b, c) @ (H(d) @ H(e))) . (omega2{0}(a, (b @ c), d) @ H(e)) . (chi1{H}((a @ b @ c @ d), e))] &
# @ [(omega2{0}(a, b, c) @ (H(d) @ H(e))) . (chi1{H}((a @ b @ c), d) @ H(e)) . (chi1{H}((a @ b @ c @ d), e))]]
omega4Target2 = ensureEqAEMol2(
    comp2s(
        comp1s(
            comp0s(_Ha, _Hb, omega2.fprimf([_H], _c, _d, _e)),
            comp0(_Ha, chi1.fprimf([_H], _b, comp0s(_c, _d, _e))),
            chi1.fprimf([_H], _a, comp0s(_b, _c, _d, _e))),
        comp1s(
            comp0s(_Ha, _Hb, chi1.fprimf([_H], _c, _d), _He),
            comp0(_Ha, omega2.fprimf([_H], _b, comp0(_c, _d), _e)),
            chi1.fprimf([_H], _a, comp0s(_b, _c, _d, _e))),
        comp1s(
            comp0s(_Ha, omega2.fprimf([_H], _b, _c, _d), _He),
            comp0(_Ha, chi1.fprimf([_H], comp0s(_b, _c, _d), _e)),
            chi1.fprimf([_H], _a, comp0s(_b, _c, _d, _e))),
        comp1s(
            comp0s(_Ha, chi1.fprimf([_H], _b, _c), _Hd, _He),
            comp0s(_Ha, chi1.fprimf([_H], comp0(_b, _c), _d), _He),
            omega2.fprimf([_H], _a, comp0s(_b, _c, _d), _e)),
        comp1s(
            comp0s(_Ha, chi1.fprimf([_H], _b, _c), _Hd, _He),
            comp0(omega2.fprimf([_H], _a, comp0(_b, _c), _d), _He),
            chi1.fprimf([_H], comp0s(_a, _b, _c, _d), _e)),
        comp1s(
            comp0s(omega2.fprimf([_H], _a, _b, _c), _Hd, _He),
            comp0(chi1.fprimf([_H], comp0s(_a, _b, _c), _d), _He),
            chi1.fprimf([_H], comp0s(_a, _b, _c, _d), _e))
    )
)


# NOTE: This has two components, giving an axiom.
omega4Paths3 = findPaths3(omega4Source2, omega4Target2, [omega3])


tensorShuf30Source0 = ensureMol0(comp0s(_Ha, _Hb, _Hc))
tensorShuf30Target0 = ensureMol0(app(_H, comp0s(_b, _c, _a)))

def _explorePaths1():
    tensorShuf30Paths1 = findPaths1(tensorShuf30Source0, tensorShuf30Target0, [shuffle1, chi1])
    exploreFromToAdjCells2(tensorShuf30Paths1, tensorShuf30Paths1, [], [omega2, u2, shuffle20, shuffle1Dim2], [omega2, u2])

#_explorePaths1()

# Too many paths; we try limiting number of shuffle1s to <= 2.  Also limit to performing all chi1s, then all shuffle1s (or vice-versa)
def numberShuffles(cell1):
    famList = lmap(lambda x: x.family, basePrimList(cell1))
    return len(list(filter(lambda x: x == shuffle1, famList)))

def basePrimList(cell1):
    l = []
    x = ensureEqMol1(cell1)
    for atom in next(iter(x.mol1s)).atom1s:
        if isinstance(atom.p1, FunctorPrim1):
            l.extend(basePrimList(atom.p1.eqMol1))
        else:
            l.append(atom.p1)
    return l

def _chi1Shuffle1NoMix(cell1):
    primList = basePrimList(cell1)
    famList = lmap(lambda x: x.family, primList)
    start = famList[0]
    switched = False
    for f in famList:
        if f != start:
            if switched:
                return False
            else:
                switched = True
            start =f
    return True

# TODO: Have to be careful about comp0(_Ha, ...) versus comp0(H(ensureEqMol1(_a), ..)
# tensorShuf30Paths1 = findPaths(tensorShuf30Source0,
#     lambda x: (x["current"] == tensorShuf30Target0) and (numberShuffles(collate1(x["cells"])) <= 2 and _chi1Shuffle1NoMix(collate1(x["cells"])) if len(x["cells"]) > 0 else False),
#     lambda x: (numberShuffles(collate1(x["cells"])) > 2) or not(_chi1Shuffle1NoMix(collate1(x["cells"]))) if len(x["cells"]) > 0 else False,
#     [shuffle1, chi1], cellsAway1, collate1)
# searchForPathPairs2AdjFam(tensorShuf30Paths1, [omega2, u2.adj, shuffle20, shuffle21, shuffle1Dim2, chi1Dim2], [])
# Take choice from Schoomer-Pries.
# [(shuffle1(H(a), (H(b) @ H(c)))) . (H(b) @ chi1{H}(c, a)) . (chi1{H}(b, (c @ a)))]
# --
# [(chi1{H}(a, b) @ H(c)) . (chi1{H}((a @ b), c)) . (H([(shuffle1(a, b) @ c)])) . (H([(b @ shuffle1(a, c))]))]
# tensorShuf30Source1 = ensureEqMol1(
#     comp1s(
#         shuffle1.fprim(_Ha, comp0(_Hb, _Hc)),
#         comp0(_Hb, chi1.fprimf([_H], _c, _a)),
#         chi1.fprimf([_H], _b, comp0(_c, _a))
#     )
# )
# tensorShuf30Target1 = ensureEqMol1(
#     comp1s(
#         comp0(chi1.fprimf([_H], _a, _b), _Hc),
#         chi1.fprimf([_H], comp0(_a, _b), _c),
#         app(_H, comp1(comp0(shuffle1.fprim(_a, _b), _c), comp0(_b, shuffle1.fprim(_a, _c))))
#     )
# )
#tensorShuf30Paths2 = findPaths2(tensorShuf30Source1, tensorShuf30Target1, [omega2, omega2.adj, u2, u2.adj, chi1Dim2, chi1Dim2.adj, shuffle20, shuffle20.adj, shuffle21, shuffle21.adj, shuffle1Dim2, shuffle1Dim2.adj])
# [[(shuffle1(H(a), (H(b) @ H(c)))) . omega2{H}(b, c, a)] &
# [shuffle1_dim2_adj([1_{H(a)}], [(chi1{H}(b, c))]) . (chi1{H}((b @ c), a))] &
# [(H(a) @ chi1{H}(b, c)) . u2_adj{H}(a, (b @ c))] &
# [omega2{H}(a, b, c) . (H([(shuffle1(a, (b @ c)))]))] &
# [(chi1{H}(a, b) @ H(c)) . (chi1{H}((a @ b), c)) . H([[shuffle20(a, b, c)]])]]
tensorShuf30Source2 = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1(
            shuffle1.fprim(_Ha, comp0(_Hb, _Hc)),
            omega2.fprimf([_H], _b, _c, _a)),
        comp1(
            shuffle1Dim2.adj.fprim(ensureEqMol1(_Ha), chi1.fprimf([_H], _b, _c)),
            chi1.fprimf([_H], comp0(_b, _c), _a)),
        comp1(
            comp0(_Ha, chi1.fprimf([_H], _b, _c)),
            u2.fprimf([_H], _a, comp0(_b, _c))),
        comp1(
            omega2.fprimf([_H], _a, _b, _c),
            app(_H, shuffle1.fprim(_a, comp0(_b, _c)))),
        comp1s(
            comp0(chi1.fprimf([_H], _a, _b), _Hc),
            chi1.fprimf([_H], comp0(_a, _b), _c),
            app(_H, shuffle20.fprim(_a, _b, _c)))
    )
), [_a, _b, _c], [_H])
# [[shuffle20(H(a), H(b), H(c)) . (H(b) @ chi1{H}(c, a)) . (chi1{H}(b, (c @ a)))] &
# [(shuffle1(H(a), H(b)) @ H(c)) . (H(b) @ shuffle1(H(a), H(c))) . omega2{H}(b, c, a)] &
# [shuffle20_adj(H(a), H(b), H(c)) . (chi1{H}(b, c) @ H(a)) . (chi1{H}((b @ c), a))] &
# [shuffle1_dim2_adj([1_{H(a)}], [(chi1{H}(b, c))]) . (chi1{H}((b @ c), a))] &
# [(H(a) @ chi1{H}(b, c)) . u2_adj{H}(a, (b @ c))] &
# @[omega2{H}(a, b, c) . (H([(shuffle1(a, (b @ c)))]))] &
# [(chi1{H}(a, b) @ H(c)) . (chi1{H}((a @ b), c)) . H([[shuffle20(a, b, c)]])]]
tensorShuf30Target2 = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1s(
            shuffle20.fprim(_Ha, _Hb, _Hc),
            comp0(_Hb, chi1.fprimf([_H], _c, _a)),
            chi1.fprimf([_H], _b, comp0(_c, _a))),
        comp1s(
            comp0(shuffle1.fprim(_Ha, _Hb), _Hc),
            comp0(_Hb, shuffle1.fprim(_Ha, _Hc)),
            omega2.fprimf([_H], _b, _c, _a)),
        comp1s(
            shuffle20.adj.fprim(_Ha, _Hb, _Hc),
            comp0(chi1.fprimf([_H], _b, _c), _Ha),
            chi1.fprimf([_H], comp0(_b, _c), _a)),
        comp1s(
            shuffle1Dim2.adj.fprim(ensureEqMol1(_Ha), chi1.fprimf([_H], _b, _c)),
            chi1.fprimf([_H], comp0(_b, _c), _a)),
        comp1(
            comp0(_Ha, chi1.fprimf([_H], _b, _c)),
            u2.fprimf([_H], _a, comp0(_b, _c))),
        comp1(
            omega2.fprimf([_H], _a, _b, _c),
            app(_H, shuffle1.fprim(_a, comp0(_b, _c)))),
        comp1s(
            comp0(chi1.fprimf([_H], _a, _b), _Hc),
            chi1.fprimf([_H], comp0(_a, _b), _c),
            app(_H, shuffle20.fprim(_a, _b, _c)))
    )
), [_a, _b, _c], [_H])
tensorShuf30 = PrimitiveFamily("tensorShuf30", 2, 1, [0, 0, 0], tensorShuf30Source2, tensorShuf30Target2)


tensorShuf31Source0 = ensureMol0(comp0s(_Ha, _Hb, _Hc))
tensorShuf31Target0 = ensureMol0(app(_H, comp0s(_b, _c, _a)))

tensorShuf31Source1 = ensureEqMol1(
    comp1s(
        shuffle1.fprim(comp0(_Ha, _Hb), _Hc),
        comp0(chi1.fprimf([_H], _c, _a), _Hb),
        chi1.fprimf([_H], comp0(_c, _a), _b)
    )
)
tensorShuf31Target1 = ensureEqMol1(
    comp1s(
        comp0(_Ha, chi1.fprimf([_H], _b, _c)),
        chi1.fprimf([_H], _a, comp0(_b, _c)),
        app(_H, comp0(_a, shuffle1.fprim(_b, _c))),
        app(_H, comp0(shuffle1.fprim(_a, _c), _b))
    )
)
#tensorShuf31Paths2 = findPaths2(tensorShuf31Source1, tensorShuf31Target1, [omega2, omega2.adj, u2, u2.adj, chi1Dim2, chi1Dim2.adj, shuffle20, shuffle20.adj, shuffle21, shuffle21.adj, shuffle1Dim2, shuffle1Dim2.adj])
ww = [omega2, omega2.adj, u2, u2.adj, chi1Dim2, chi1Dim2.adj, shuffle20, shuffle20.adj, shuffle21, shuffle21.adj, shuffle1Dim2, shuffle1Dim2.adj]


# We pick id_{H(x \pt y)} instead of H(id_{x \pt y}) as part of target for matching ease reasons.
u3Source1 = ensureEqMol1(
    comp1s(
        shuffle1.fprim(_Ha, _Hb),
        shuffle1.fprim(_Hb, _Ha),
        chi1.fprimf([_H], _a, _b)
    )
)
u3Target1 = ensureEqMol1(
    comp1s(
        chi1.fprimf([_H], _a, _b)
    )
)
#u3Paths2 = findPaths2(u3Source1, u3Target1, [syllepsis2, u2, appComp1, appIdentity1.adj])
# [[syllepsis2(H(a), H(b)) . (chi1{H}(a, b))]]
u3Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp1(
        syllepsis2.fprim(_Ha, _Hb),
        chi1.fprimf([_H], _a, _b)
    )
), [_a, _b], [_H])
# [[(shuffle1(H(a), H(b))) . u2{H}(b, a)] &
# [u2{H}(a, b) . (H([(shuffle1(b, a))]))] &
# [(chi1{H}(a, b)) . appComp1{H}([(shuffle1(a, b))], [(shuffle1(b, a))])] &
# [(chi1{H}(a, b)) . H([[syllepsis2(a, b)]])] &
# [(chi1{H}(a, b)) . appIdentity1_adj{H}((a @ b))]]
u3Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1s(
            shuffle1.fprim(_Ha, _Hb),
            u2.fprimf([_H], _b, _a)),
        comp1s(
            u2.fprimf([_H], _a, _b),
            app(_H, shuffle1.fprim(_b, _a))),
        comp1s(
            chi1.fprimf([_H], _a, _b),
            appComp1.fprimf([_H], shuffle1.fprim(_a, _b), shuffle1.fprim(_b, _a))),
        comp1s(
            chi1.fprimf([_H], _a, _b),
            app(_H, syllepsis2.fprim(_a, _b))),
        comp1s(
            chi1.fprimf([_H], _a, _b),
            appIdentity1.adj.fprimf([_H], comp0(_a, _b)))
    )
), [_a, _b], [_H])
u3 = PrimitiveFamily("u3", 3, 1, [0, 0], u3Source, u3Target)
