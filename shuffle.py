from fill import *
from ontology import *
from primfamily import *
from transfor import *

# Take two lists, and returns a list of all the shuffles of one into the other.
# Gives them in order such that inserting all at the beginning is the first element, and inserting all at the end is the last.
def shuffles(a, b):
    if len(a) == 0:
        return [b]
    if len(b) == 0:
        return [a]
    return list(map(lambda z: [a[0]] + z, shuffles(a[1:], b))) + list(map(lambda z: [b[0]] + z, shuffles(a, b[1:])))

_a = ConstPrim0("a")
_b = ConstPrim0("b")
_c = ConstPrim0("c")
_d = ConstPrim0("d")

shuffle1Source = minimalASTFromMol0(comp0(_a, _b), [_a, _b])
shuffle1Target = minimalASTFromMol0(comp0(_b, _a), [_a, _b])
shuffle1 = PrimitiveFamily("shuffle1", 1, 0, [0, 0], shuffle1Source, shuffle1Target)
shuffle1Dim2 = tritrans2CellPrimFamily(shuffle1)
shuffle1Dim3 = tritrans3CellPrimFamily(shuffle1, shuffle1Dim2)
shuffle1Pi = tritransPiPrimFamily(shuffle1, shuffle1Dim2)
shuffle1Pi.adj.isDegen = lambda params: (isIdMol1(params[0]) and isIdMol1(params[1])) or (isIdMol1(params[2]) and isIdMol1(params[3]))
shuffle1M = tritransMPrimFamily(shuffle1, shuffle1Dim2)

# We want paths involving only shuffled vertices, so we need to use a custom endFunc
#shuffle20Vertices = list(map(lambda x: comp0s(*x), shuffles([_a], [_b, _c])))
#shuffle20Paths1 = findPaths(shuffle20Vertices[0], targetEndCond(shuffle20Vertices[-1]), lambda x: x["current"] not in shuffle20Vertices, [shuffle1], cellsAway1, collate1)
# [(shuffle1(a, (b @ c)))] -- [(shuffle1(a, b) @ c) . (b @ shuffle1(a, c))]
shuffle20Source = minimalASTFromEqMol1(ensureEqMol1(
    shuffle1.fprim(_a, comp0(_b, _c))), [_a, _b, _c])
shuffle20Target = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(shuffle1.fprim(_a, _b), _c),
        comp0(_b, shuffle1.fprim(_a, _c)))), [_a, _b, _c])
shuffle20 = PrimitiveFamily("shuffle20", 2, 0, [0, 0, 0], shuffle20Source, shuffle20Target)
shuffle20Dim3 = trimod3CellPrimFamily(shuffle20, {shuffle1: shuffle1Dim2})

#shuffle21Vertices = list(map(lambda x: comp0s(*x), shuffles([_a, _b], [_c])))
#shuffle21Paths1 = findPaths(shuffle21Vertices[0], targetEndCond(shuffle21Vertices[-1]), lambda x: x["current"] not in shuffle21Vertices, [shuffle1], cellsAway1, collate1)
# [(shuffle1((a @ b), c))] -- [(a @ shuffle1(b, c)) . (shuffle1(a, c) @ b)]
shuffle21Source = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(_a, shuffle1.fprim(_b, _c)),
        comp0(shuffle1.fprim(_a, _c), _b))), [_a, _b, _c])
shuffle21Target = minimalASTFromEqMol1(ensureEqMol1(
    shuffle1.fprim(comp0(_a, _b), _c)), [_a, _b, _c])
shuffle21 = PrimitiveFamily("shuffle21", 2, 0, [0, 0, 0], shuffle21Source, shuffle21Target)
shuffle21Dim3 = trimod3CellPrimFamily(shuffle21, {shuffle1: shuffle1Dim2})

# NOTE: For Dim1 in all the below, We picked the codomain candidate that had no tensors inside of the individual shuffles, and the domain candidate that was a single prim.
# ABCD -> BCDA (n = 1, k = 3, insert A into BCD)
#shuffle30Vertices = list(map(lambda x: comp0s(*x), shuffles([_a], [_b, _c, _d])))
#shuffle30Paths1 = findPaths(shuffle30Vertices[0], targetEndCond(shuffle30Vertices[-1]), lambda x: x["current"] not in shuffle30Vertices, [shuffle1], cellsAway1, collate1)
#searchForPathPairs2(shuffle30Paths1, [shuffle20, shuffle21])
# [(shuffle1(a, (b @ c @ d)))] -- [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, c) @ d) . ((b @ c) @ shuffle1(a, d))]
#shuffle30Source1 = ensureEqMol1(shuffle1.fprim(_a, comp0s(_b, _c, _d)))
#shuffle30Target1 = ensureEqMol1(comp1s(
#     comp0s(shuffle1.fprim(_a, _b), _c, _d),
#     comp0s(_b, shuffle1.fprim(_a, _c), _d),
#     comp0s(_b, _c, shuffle1.fprim(_a, _d))
#))
# shuffle30Paths2 = findPaths2(shuffle30Source1, shuffle30Target1, [shuffle20])
# [[shuffle20(a, (b @ c), d)] & [(shuffle20(a, b, c) @ d) . ((b @ c) @ shuffle1(a, d))]]  --  [[shuffle20(a, b, (c @ d))] & [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle20(a, c, d))]]
shuffle30Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        shuffle20.fprim(_a, comp0(_b, _c), _d),
        comp1(
            comp0(shuffle20.fprim(_a, _b, _c), _d),
            comp0s(_b, _c, shuffle1.fprim(_a, _d))
        )
    )), [_a, _b, _c, _d])
shuffle30Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        shuffle20.fprim(_a, _b, comp0(_c, _d)),
        comp1(
            comp0s(shuffle1.fprim(_a, _b), _c, _d),
            comp0(_b, shuffle20.fprim(_a, _c, _d))
        )
    )), [_a, _b, _c, _d])
shuffle30 = PrimitiveFamily("shuffle30", 3, 0, [0, 0, 0, 0], shuffle30Source, shuffle30Target)

# ABCD -> CDAB (n = 2, k = 2, insert AB into CD)
#shuffle31Vertices = list(map(lambda x: comp0s(*x), shuffles([_a, _b], [_c, _d])))
#shuffle31Paths1 = findPaths(shuffle31Vertices[0], targetEndCond(shuffle31Vertices[-1]), lambda x: x["current"] not in shuffle31Vertices, [shuffle1], cellsAway1, collate1)
#searchForPathPairs2(shuffle31Paths1, [shuffle20, shuffle21])
#[(shuffle1((a @ b), (c @ d)))] -- [(a @ shuffle1(b, c) @ d) . (shuffle1(a, c) @ (b @ d)) . ((c @ a) @ shuffle1(b, d)) . (c @ shuffle1(a, d) @ b)]
#shuffle31Source1  = ensureEqMol1(shuffle1.fprim(comp0(_a, _b), comp0(_c, _d)))
#shuffle31Target1 = ensureEqMol1(comp1s(
#     comp0s(_a, shuffle1.fprim(_b, _c), _d),
#     comp0s(shuffle1.fprim(_a,_c), _b, _d),
#     comp0s(_c, _a, shuffle1.fprim(_b, _d)),
#     comp0s(_c, shuffle1.fprim(_a, _d), _b)
# ))
#shuffle31Paths2 = findPaths2(shuffle31Source1, shuffle31Target1, [shuffle20, shuffle21])
# [[shuffle20((a @ b), c, d)] & [(shuffle21(a, b, c) @ d) . (c @ shuffle1((a @ b), d))] & [(a @ shuffle1(b, c) @ d) . (shuffle1(a, c) @ (b @ d)) . (c @ shuffle21(a, b, d))]]
# -- [[shuffle21(a, b, (c @ d))] & [(a @ shuffle20(b, c, d)) . (shuffle1(a, (c @ d)) @ b)] & [(a @ shuffle1(b, c) @ d) . ((a @ c) @ shuffle1(b, d)) . (shuffle20(a, c, d) @ b)]]
# shuffle31Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         shuffle20.fprim(comp0(_a, _b), _c, _d),
#         comp1(
#             comp0(shuffle21.fprim(_a, _b, _c), _d),
#             comp0(_c, shuffle1.fprim(comp0(_a, _b), _d))
#         ),
#         comp1s(
#             comp0s(_a, shuffle1.fprim(_b, _c), _d),
#             comp0s(shuffle1.fprim(_a,_c), _b, _d),
#             comp0(_c, shuffle21.fprim(_a, _b, _d))
#         )
#     )), [_a, _b, _c, _d])
# shuffle31Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         shuffle21.fprim(_a, _b, comp0(_c, _d)),
#         comp1(
#             comp0(_a, shuffle20.fprim(_b, _c, _d)),
#             comp0(shuffle1.fprim(_a, comp0(_c, _d)), _b)
#         ),
#         comp1s(
#             comp0s(_a, shuffle1.fprim(_b, _c), _d),
#             comp0s(_a, _c, shuffle1.fprim(_b, _d)),
#             comp0(shuffle20.fprim(_a, _c, _d), _b)
#         )
#     )), [_a, _b, _c, _d])
# shuffle31 = PrimitiveFamily("shuffle31", 3, 0, [0, 0, 0, 0], shuffle31Source, shuffle31Target)

# ABCD -> DABC (n = 3, k = 1, insert ABC into D)
# shuffle32Vertices = list(map(lambda x: comp0s(*x), shuffles([_a, _b, _c], [_d])))
# shuffle32Paths1 = findPaths(shuffle32Vertices[0], targetEndCond(shuffle32Vertices[-1]), lambda x: x["current"] not in shuffle32Vertices, [shuffle1], cellsAway1, collate1)
# searchForPathPairs2(shuffle32Paths1, [shuffle20, shuffle21])
# [(shuffle1((a @ b @ c), d))] -- [((a @ b) @ shuffle1(c, d)) . (a @ shuffle1(b, d) @ c) . (shuffle1(a, d) @ (b @ c))]
# shuffle32Source1 = ensureEqMol1(shuffle1.fprim(comp0s(_a, _b, _c), _d))
# shuffle32Target1 = ensureEqMol1(
#     comp1s(
#         comp0s(_a, _b, shuffle1.fprim(_c, _d)),
#         comp0s(_a, shuffle1.fprim(_b, _d), _c),
#         comp0s(shuffle1.fprim(_a, _d), _b, _c)
#     ))
#shuffle32Paths2 = findPaths2(shuffle32Source1, shuffle32Target1, [shuffle20, shuffle21])
# [[shuffle21((a @ b), c, d)] & [((a @ b) @ shuffle1(c, d)) . (shuffle21(a, b, d) @ c)]]
# -- [[shuffle21(a, (b @ c), d)] & [(a @ shuffle21(b, c, d)) . (shuffle1(a, d) @ (b @ c))]]
# shuffle32Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         shuffle21.fprim(comp0(_a, _b), _c, _d),
#         comp1(
#             comp0s(_a, _b, shuffle1.fprim(_c, _d)),
#             comp0(shuffle21.fprim(_a, _b, _d), _c)
#         )
#     )), [_a, _b, _c, _d])
# shuffle32Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         shuffle21.fprim(_a, comp0(_b, _c), _d),
#         comp1(
#             comp0(_a, shuffle21.fprim(_b, _c, _d)),
#             comp0s(shuffle1.fprim(_a, _d), _b, _c)
#         )
#     )), [_a, _b, _c ,_d])
# shuffle32 = PrimitiveFamily("shuffle32", 3, 0, [0, 0, 0, 0], shuffle32Source, shuffle32Target)


#breen3Paths1 = findPaths1(comp0s(_a, _b, _c), comp0s(_c, _b, _a), [shuffle1])
# This choice is different from Schommer-Pries, but has the same number of cells , and the 2-cell paths involve fewer adjoints.
# [(a @ shuffle1(b, c)) . (shuffle1(a, (c @ b)))] -- [(shuffle1((a @ b), c)) . (c @ shuffle1(a, b))]
# breen3Source1 = ensureEqMol1(comp1s(
#     comp0(_a, shuffle1.fprim(_b, _c)),
#     shuffle1.fprim(_a, comp0(_c, _b))))
# breen3Target1 = ensureEqMol1(comp1s(
#     shuffle1.fprim(comp0(_a, _b), _c),
#     comp0(_c, shuffle1.fprim(_a, _b))))
#breen3Paths2 = findPaths2(breen3Source1, breen3Target1, [shuffle20,  shuffle21.adj, shuffle1Dim2])
# [[(a @ shuffle1(b, c)) . shuffle20(a, c, b)] &
# [shuffle21_adj(a, b, c) . (c @ shuffle1(a, b))]]
# -- [[shuffle1_dim2([1_{a}], [(shuffle1(b, c))])] &
# [shuffle20(a, b, c) . (shuffle1(b, c) @ a)] &
# [(shuffle1(a, b) @ c) . shuffle21_adj(b, a, c)] &
# [shuffle1_dim2([(shuffle1(a, b))], [1_{c}])]]
# breen3Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         comp1s(
#             comp0(_a, shuffle1.fprim(_b, _c)),
#             shuffle20.fprim(_a, _c, _b)),
#         comp1s(
#             shuffle21.adj.fprim(_a, _b, _c),
#             comp0(_c, shuffle1.fprim(_a, _b))))), [_a, _b, _c])
# breen3Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         shuffle1Dim2.fprim(ensureEqMol1(_a), shuffle1.fprim(_b, _c)),
#         comp1s(
#             shuffle20.fprim(_a, _b, _c),
#             comp0(shuffle1.fprim(_b, _c), _a)),
#         comp1s(
#             comp0(shuffle1.fprim(_a, _b), _c),
#             shuffle21.adj.fprim(_b, _a, _c)),
#         shuffle1Dim2.fprim(shuffle1.fprim(_a, _b), ensureEqMol1(_c)))), [_a, _b, _c])
# breen3 = PrimitiveFamily("breen3", 3, 0, [0, 0, 0], breen3Source, breen3Target)


# Alternative
# # [(shuffle1(a, (b @ c))) . (b @ shuffle1(c, a)) . (shuffle1(b, (a @ c))) . (shuffle1(a, (c @ b)))]
# # -- [(shuffle1(a, b) @ c) . (b @ shuffle1(a, c)) . (shuffle1(b, (c @ a))) . (shuffle1(c, a) @ b) . (shuffle1(a, (c @ b)))]
# breen3Source1 = ensureEqMol1(comp1s(
#     shuffle1.fprim(_a, comp0(_b, _c)),
#     comp0(_b, shuffle1.fprim(_c, _a)),
#     shuffle1.fprim(_b, comp0(_a, _c)),
#     shuffle1.fprim(_a, comp0(_c, _b))
# ))
# breen3Target1 = ensureEqMol1(comp1s(
#     comp0(shuffle1.fprim(_a, _b), _c),
#     comp0(_b, shuffle1.fprim(_a, _c)),
#     shuffle1.fprim(_b, comp0(_c, _a)),
#     comp0(shuffle1.fprim(_c, _a), _b),
#     shuffle1.fprim(_a, comp0(_c, _b))
# ))
# breen3Paths2 = findPaths2(breen3Source1, breen3Target1, [shuffle20, shuffle21, shuffle1Dim2, shuffle1Dim2.adj])
# # [[shuffle20(a, b, c) . (b @ shuffle1(c, a)) . (shuffle1(b, (a @ c))) . (shuffle1(a, (c @ b)))] &
# # [(shuffle1(a, b) @ c) . (b @ shuffle1(a, c)) . shuffle1_dim2([1_{b}], [(shuffle1(c, a))]) . (shuffle1(a, (c @ b)))]]
# # --
# # [[shuffle20(a, b, c) . (b @ shuffle1(c, a)) . (shuffle1(b, (a @ c))) . (shuffle1(a, (c @ b)))] &
# # [(shuffle1(a, b) @ c) . shuffle1_dim2([1_{b}], [(shuffle1(a, c)) . (shuffle1(c, a))]) . (shuffle1(a, (c @ b)))] &
# # [(shuffle1(a, b) @ c) . shuffle1_dim2_adj([1_{b}], [(shuffle1(a, c))]) . (shuffle1(c, a) @ b) . (shuffle1(a, (c @ b)))]]
# breen3Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         comp1s(
#             shuffle20.fprim(_a, _b, _c),
#             comp0(_b, shuffle1.fprim(_c, _a)),
#             shuffle1.fprim(_b, comp0(_a, _c)),
#             shuffle1.fprim(_a, comp0(_c, _b))),
#         comp1s(
#             comp0(shuffle1.fprim(_a, _b), _c),
#             comp0(_b, shuffle1.fprim(_a, _c)),
#             shuffle1Dim2.fprim(ensureEqMol1(_b), ensureEqMol1(shuffle1.fprim(_c, _a))),
#             shuffle1.fprim(_a, comp0(_c, _b))))), [_a, _b, _c])
# breen3Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
#     comp2s(
#         comp1s(
#             shuffle20.fprim(_a, _b, _c),
#             comp0(_b, shuffle1.fprim(_c, _a)),
#             shuffle1.fprim(_b, comp0(_a, _c)),
#             shuffle1.fprim(_a, comp0(_c, _b))),
#         comp1s(
#             comp0(shuffle1.fprim(_a, _b), _c),
#             shuffle1Dim2.fprim(ensureEqMol1(_b), comp1s(shuffle1.fprim(_a, _c), shuffle1.fprim(_c, _a))),
#             shuffle1.fprim(_a, comp0(_c, _b))),
#         comp1s(
#             comp0(shuffle1.fprim(_a, _b), _c),
#             shuffle1Dim2.adj.fprim(ensureEqMol1(_b), ensureEqMol1(shuffle1.fprim(_a, _c))),
#             comp0(shuffle1.fprim(_c, _a), _b),
#             shuffle1.fprim(_a, comp0(_c, _b))))), [_a, _b, _c])

# ALTERNATIVE EXPLORATION
# Since there is an expansion of breen3 using only shuffle20 and shuffle21, redo shuffle31 and shuffle32 to also use only these.
# This should give us best expansion of breen4.


# Shuffle30 is same

# ABCD -> CDAB (n = 2, k = 2, insert AB into CD)
#shuffle31Vertices = list(map(lambda x: comp0s(*x), shuffles([_a, _b], [_c, _d])))
#shuffle31Paths1 = findPaths(shuffle31Vertices[0], targetEndCond(shuffle31Vertices[-1]), lambda x: x["current"] not in shuffle31Vertices, [shuffle1], cellsAway1, collate1)
#searchForPathPairs2(shuffle31Paths1, [shuffle20, shuffle21])
#[(a @ shuffle1(b, (c @ d))) . (shuffle1(a, (c @ d)) @ b)] -- [(shuffle1((a @ b), c) @ d) . (c @ shuffle1((a @ b), d))]
# shuffle31Source1  = ensureEqMol1(
#     comp1(
#         comp0(_a, shuffle1.fprim(_b, comp0(_c, _d))),
#         comp0(shuffle1.fprim(_a, comp0(_c, _d)), _b)))
# shuffle31Target1 = ensureEqMol1(
#     comp1(
#         comp0(shuffle1.fprim(comp0(_a, _b), _c), _d),
#         comp0(_c, shuffle1.fprim(comp0(_a, _b), _d))))
# shuffle31Paths2 = findPaths2(shuffle31Source1, shuffle31Target1, [shuffle20, shuffle21])
# # [[shuffle21(a, b, (c @ d))] &
# # [shuffle20((a @ b), c, d)]]
# # -- [[(a @ shuffle20(b, c, d)) . (shuffle1(a, (c @ d)) @ b)] &
# # [(a @ shuffle1(b, c) @ d) . ((a @ c) @ shuffle1(b, d)) . (shuffle20(a, c, d) @ b)] &
# # [(shuffle21(a, b, c) @ d) . ((c @ a) @ shuffle1(b, d)) . (c @ shuffle1(a, d) @ b)] &
# # [(shuffle1((a @ b), c) @ d) . (c @ shuffle21(a, b, d))]]
shuffle31Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        shuffle21.fprim(_a, _b, comp0(_c, _d)),
        shuffle20.fprim(comp0(_a, _b), _c, _d)
    )), [_a, _b, _c, _d])
shuffle31Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1s(
            comp0(_a, shuffle20.fprim(_b, _c, _d)),
            comp0(shuffle1.fprim(_a, comp0(_c, _d)), _b)),
        comp1s(
            comp0s(_a, shuffle1.fprim(_b, _c), _d),
            comp0s(_a, _c, shuffle1.fprim(_b, _d)),
            comp0(shuffle20.fprim(_a, _c, _d), _b)),
        comp1s(
            comp0(shuffle21.fprim(_a, _b, _c), _d),
            comp0s(_c, _a, shuffle1.fprim(_b, _d)),
            comp0s(_c, shuffle1.fprim(_a, _d), _b)),
        comp1s(
            comp0(shuffle1.fprim(comp0(_a, _b), _c), _d),
            comp0(_c, shuffle21.fprim(_a, _b, _d)))
    )), [_a, _b, _c, _d])
shuffle31 = PrimitiveFamily("shuffle31", 3, 0, [0, 0, 0, 0], shuffle31Source, shuffle31Target)
#
# # ABCD -> DABC (n = 3, k = 1, insert ABC into D)
# # shuffle32Vertices = list(map(lambda x: comp0s(*x), shuffles([_a, _b, _c], [_d])))
# # shuffle32Paths1 = findPaths(shuffle32Vertices[0], targetEndCond(shuffle32Vertices[-1]), lambda x: x["current"] not in shuffle32Vertices, [shuffle1], cellsAway1, collate1)
# # searchForPathPairs2(shuffle32Paths1, [shuffle20, shuffle21])
# # [((a @ b) @ shuffle1(c, d)) . (a @ shuffle1(b, d) @ c) . (shuffle1(a, d) @ (b @ c))] -- [(shuffle1((a @ b @ c), d))]
# shuffle32Source1 = ensureEqMol1(
#     comp1s(
#         comp0s(_a, _b, shuffle1.fprim(_c, _d)),
#         comp0s(_a, shuffle1.fprim(_b, _d), _c),
#         comp0s(shuffle1.fprim(_a, _d), _b, _c)
#     ))
# shuffle32Target1 = ensureEqMol1(shuffle1.fprim(comp0s(_a, _b, _c), _d))
# shuffle32Paths2 = findPaths2(shuffle32Source1, shuffle32Target1, [shuffle20, shuffle21j])
# # [[shuffle21((a @ b), c, d)] & [((a @ b) @ shuffle1(c, d)) . (shuffle21(a, b, d) @ c)]]
# # -- [[shuffle21(a, (b @ c), d)] & [(a @ shuffle21(b, c, d)) . (shuffle1(a, d) @ (b @ c))]]
shuffle32Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1(
            comp0(_a, shuffle21.fprim(_b, _c, _d)),
            comp0s(shuffle1.fprim(_a, _d), _b, _c)
        ),
        shuffle21.fprim(_a, comp0(_b, _c), _d),
    )), [_a, _b, _c ,_d])
shuffle32Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1(
            comp0s(_a, _b, shuffle1.fprim(_c, _d)),
            comp0(shuffle21.fprim(_a, _b, _d), _c)
        ),
        shuffle21.fprim(comp0(_a, _b), _c, _d)
    )), [_a, _b, _c, _d])
shuffle32 = PrimitiveFamily("shuffle32", 3, 0, [0, 0, 0, 0], shuffle32Source, shuffle32Target)
#
#
# #breen3Paths1 = findPaths1(comp0s(_a, _b, _c), comp0s(_c, _b, _a), [shuffle1])
# # This choice is different from Schommer-Pries, but has the same number of cells , and the 2-cell paths involve fewer adjoints.
# # [(a @ shuffle1(b, c)) . (shuffle1(a, (c @ b)))] -- [(shuffle1((a @ b), c)) . (c @ shuffle1(a, b))]
# # breen3Source1 = ensureEqMol1(comp1s(
# #     comp0(_a, shuffle1.fprim(_b, _c)),
# #     shuffle1.fprim(_a, comp0(_c, _b))))
# # breen3Target1 = ensureEqMol1(comp1s(
# #     shuffle1.fprim(comp0(_a, _b), _c),
# #     comp0(_c, shuffle1.fprim(_a, _b))))
# # breen3Paths2 = findPaths2(breen3Source1, breen3Target1, [shuffle20,  shuffle21, shuffle1Dim2])
# # [[(a @ shuffle1(b, c)) . shuffle20(a, c, b)] &
# # [shuffle21(a, b, c) . (c @ shuffle1(a, b))]]
# # -- [[shuffle1_dim2([1_{a}], [(shuffle1(b, c))])] &
# # [shuffle20(a, b, c) . (shuffle1(b, c) @ a)] &
# # [(shuffle1(a, b) @ c) . shuffle21(b, a, c)] &
# # [shuffle1_dim2([(shuffle1(a, b))], [1_{c}])]]
breen3Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        comp1s(
            comp0(_a, shuffle1.fprim(_b, _c)),
            shuffle20.fprim(_a, _c, _b)),
        comp1s(
            shuffle21.fprim(_a, _b, _c),
            comp0(_c, shuffle1.fprim(_a, _b))))), [_a, _b, _c])
breen3Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s(
        shuffle1Dim2.fprim(ensureEqMol1(_a), shuffle1.fprim(_b, _c)),
        comp1s(
            shuffle20.fprim(_a, _b, _c),
            comp0(shuffle1.fprim(_b, _c), _a)),
        comp1s(
            comp0(shuffle1.fprim(_a, _b), _c),
            shuffle21.fprim(_b, _a, _c)),
        shuffle1Dim2.fprim(shuffle1.fprim(_a, _b), ensureEqMol1(_c)))), [_a, _b, _c])
breen3 = PrimitiveFamily("breen3", 3, 0, [0, 0, 0], breen3Source, breen3Target)




# NOTE: Full path search here is too large, so we find a simple path and then correct it.
# Takes too long: breen4Paths1 = findPaths1(comp0s(_a, _b, _c, _d), comp0s(_d, _c, _b, _a), [shuffle1])
# def isNotSingleBraid(path):
#     if len(path["cells"]) == 0:
#         return False
#     lastCell = path["cells"][-1]
#     if len(lastCell.p1.params[0]) != 1 or len(lastCell.p1.params[1]) != 1:
#         return True
#     return False
# Too many paths to check them all; but we can check the shortest ones.
#breen4Paths1 = findPaths(comp0s(_a, _b, _c, _d), targetEndCond(comp0s(_d, _c, _b, _a)), lambda x: len(x["cells"]) > 2, [shuffle1], cellsAway1, collate1)
#searchForPathPairs2(breen4Paths1, [shuffle20,  shuffle21, shuffle1Dim2])
# Two short candidates:
# [(a @ b @ shuffle1(c, d) . (a @ shuffle1(b, (d @ c))) . (shuffle1(a, (d @ c @ b)))]
# -- [(shuffle1((a @ b @ c), d)) . (d @ shuffle1((a @ b), c)) . ((d @ c) @ shuffle1(a, b))]
breen4Source1 = fEqMol1(
    comp1s(
        comp0s(_a, _b, shuffle1.fprim(_c, _d)),
        comp0(_a, shuffle1.fprim(_b, comp0(_d, _c))),
        shuffle1.fprim(_a, comp0s(_d, _c, _b))
    )
)
breen4Target1 = fEqMol1(
    comp1s(
        shuffle1.fprim(comp0s(_a, _b, _c), _d),
        comp0(_d, shuffle1.fprim(comp0(_a, _b), _c)),
        comp0s(_d, _c, shuffle1.fprim(_a, _b))
    )
)

def _explorePaths1():
    for p1 in breen4Paths1:
        for p2 in breen4Paths1:
            paths2 = findPaths2(p1, p2, [shuffle20, shuffle21, shuffle1Dim2])
            exploreFromToAdjCells3(paths2, paths2, [shuffle30, shuffle31, shuffle32, breen3, shuffle20Dim3, shuffle21Dim3, shuffle1Dim3, shuffle1Pi], [],
                    [shuffle30, shuffle30.adj, shuffle31, shuffle31.adj, shuffle32, shuffle32.adj, breen3, breen3.adj])

#_explorePaths1()

#breen4Paths2 = list(findPaths(breen4Source1, targetEndCond(breen4Target1), lambda x: len(x["cells"]) >= 100, [shuffle20, shuffle21, shuffle1Dim2], cellsAway2, collate2))
#breen4Paths2.sort(key = len)
# breen3, shuffle31
#circleAround3(breen4Paths2[0], [shuffle30.adj, shuffle31, shuffle32.adj, breen3, shuffle20Dim3, shuffle21Dim3.adj.getMate(), shuffle1Dim3])
# exploreFromToAdjCells3(breen4Paths2, breen4Paths2, [],
#     [breen3, breen3.adj, shuffle30, shuffle31, shuffle32,  shuffle30.adj, shuffle31.adj, shuffle32.adj,
#     shuffle1Dim3, shuffle1Dim3.adj, shuffle1Dim3.getMate(), shuffle1Dim3.getMate().adj,
#     shuffle20Dim3, shuffle20Dim3.adj, shuffle20Dim3.getMate(), shuffle20Dim3.getMate().adj,
#     shuffle21Dim3, shuffle21Dim3.adj, shuffle21Dim3.getMate(), shuffle21Dim3.getMate().adj,
#     shuffle1Pi, shuffle1Pi.adj, shuffle1Pi.getMate(), shuffle1Pi.getMate().adj,
#     shuffle1Dim2.adjId1, shuffle1Dim2.adjId2, shuffle20.adjId1, shuffle20.adjId2, shuffle21.adjId1, shuffle21.adjId2], [shuffle30, shuffle31, shuffle31.adj, shuffle32.adj, breen3])
# exploreFromToAdjCells3(breen4Paths2, breen4Paths2, [breen3],
#     [shuffle30, shuffle31, shuffle31.adj, shuffle32.adj,
#     shuffle1Dim3, shuffle1Dim3.adj, shuffle1Dim3.getMate(), shuffle1Dim3.getMate().adj,
#     shuffle20Dim3, shuffle20Dim3.adj, shuffle20Dim3.getMate(), shuffle20Dim3.getMate().adj,
#     shuffle21Dim3, shuffle21Dim3.adj, shuffle21Dim3.getMate(), shuffle21Dim3.getMate().adj,
#     shuffle1Pi, shuffle1Pi.adj, shuffle1Pi.getMate(), shuffle1Pi.getMate().adj,
#     shuffle1Dim2.adjId1, shuffle1Dim2.adjId2, shuffle20.adjId1, shuffle20.adjId2, shuffle21.adjId1, shuffle21.adjId2], [shuffle30, shuffle31, shuffle31.adj, shuffle32.adj, breen3])

# ['shuffle30', 'shuffle31', 'shuffle32_adj', 'breen3_adj', 'shuffle1_dim3', 'shuffle20_dim3', 'shuffle21_dim3']

# [[shuffle1_dim2([1_{a}], [(b @ shuffle1(c, d)) . (shuffle1(b, (d @ c)))])] &
# [shuffle20(a, b, (c @ d)) . (b @ shuffle1(c, d) @ a) . (shuffle1(b, (d @ c)) @ a)] &
# [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, (c @ d))) . (shuffle1_dim2([1_{b}], [(shuffle1(c, d))]) @ a)] &
# [(shuffle1(a, b) @ (c @ d)) . shuffle21_adj(b, a, (c @ d)) . (shuffle1(c, d) @ (b @ a))] &
# [(shuffle1(a, b) @ (c @ d)) . shuffle20((b @ a), c, d) . (shuffle1(c, d) @ (b @ a))] &
# [(shuffle1_dim2([(shuffle1(a, b))], [1_{c}]) @ d) . (c @ shuffle1((b @ a), d)) . (shuffle1(c, d) @ (b @ a))] &
# [(shuffle1((a @ b), c) @ d) . (c @ shuffle1(a, b) @ d) . shuffle21_adj(c, (b @ a), d)] &
# [shuffle1_dim2([(shuffle1((a @ b), c)) . (c @ shuffle1(a, b))], [1_{d}])]]
# breen4Source2 = ensureEqAEMol2(
#     comp2s(
#         shuffle1Dim2.fprim(ensureEqMol1(_a), comp1s(comp0(_b, shuffle1.fprim(_c, _d)), shuffle1.fprim(_b, comp0(_d, _c)))),
#         comp1s(
#             shuffle20.fprim(_a, _b, comp0(_c, _d)),
#             comp0s(_b, shuffle1.fprim(_c, _d), _a),
#             comp0(shuffle1.fprim(_b, comp0(_d, _c)), _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             comp0(_b, shuffle1.fprim(_a, comp0(_c, _d))),
#             comp0(shuffle1Dim2.fprim(ensureEqMol1(_b), shuffle1.fprim(_c, _d)), _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             shuffle21.adj.fprim(_b, _a, comp0(_c, _d)),
#             comp0s(shuffle1.fprim(_c, _d), _b, _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             shuffle20.fprim(comp0(_b, _a), _c, _d),
#             comp0s(shuffle1.fprim(_c, _d), _b, _a)),
#         comp1s(
#             comp0(shuffle1Dim2.fprim(shuffle1.fprim(_a, _b), ensureEqMol1(_c)), _d),
#             comp0(_c, shuffle1.fprim(comp0(_b, _a), _d)),
#             comp0s(shuffle1.fprim(_c, _d), _b, _a)),
#         comp1s(
#             comp0(shuffle1.fprim(comp0(_a, _b), _c), _d),
#             comp0s(_c, shuffle1.fprim(_a, _b), _d),
#             shuffle21.adj.fprim(_c, comp0(_b, _a), _d)),
#         shuffle1Dim2.fprim(comp1(shuffle1.fprim(comp0(_a, _b), _c), comp0(_c, shuffle1.fprim(_a, _b))), ensureEqMol1(_d))
#     )
# )

# --
# [[shuffle1_dim2([1_{a}], [(b @ shuffle1(c, d)) . (shuffle1(b, (d @ c)))])] &
# [shuffle20(a, b, (c @ d)) . (b @ shuffle1(c, d) @ a) . (shuffle1(b, (d @ c)) @ a)] &
# @[(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, (c @ d))) . (shuffle1_dim2([1_{b}], [(shuffle1(c, d))]) @ a)] &
# [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle20(a, c, d)) . (shuffle1(b, (c @ d)) @ a) . (shuffle1(c, d) @ (b @ a))] &
# [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, c) @ d) . ((b @ c) @ shuffle1(a, d)) . (shuffle20(b, c, d) @ a) . (shuffle1(c, d) @ (b @ a))] &
# [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, c) @ d) . (shuffle1(b, c) @ (a @ d)) . (c @ shuffle21_adj(b, a, d)) . (shuffle1(c, d) @ (b @ a))] &
# [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, c) @ d) . (shuffle1(b, c) @ (a @ d)) . shuffle21_adj(c, (b @ a), d)] &
# [shuffle1_dim2([(shuffle1(a, b) @ c) . (b @ shuffle1(a, c)) . (shuffle1(b, c) @ a)], [1_{d}])] &
# [(shuffle1((a @ b @ c), d)) . (d @ shuffle1(a, b) @ c) . (d @ shuffle21_adj(b, a, c))] &
# [(shuffle1((a @ b @ c), d)) . (d @ shuffle1_dim2([(shuffle1(a, b))], [1_{c}]))]]
# breen4Target2 = ensureEqAEMol2(
#     comp2s(
#         shuffle1Dim2.fprim(ensureEqMol1(_a), comp1s(comp0(_b, shuffle1.fprim(_c, _d)), shuffle1.fprim(_b, comp0(_d, _c)))),
#         comp1s(
#             shuffle20.fprim(_a, _b, comp0(_c, _d)),
#             comp0s(_b, shuffle1.fprim(_c, _d), _a),
#             comp0(shuffle1.fprim(_b, comp0(_d, _c)), _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             comp0(_b, shuffle1.fprim(_a, comp0(_c, _d))),
#             comp0(shuffle1Dim2.fprim(ensureEqMol1(_b), shuffle1.fprim(_c, _d)), _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             comp0(_b, shuffle20.fprim(_a, _c, _d)),
#             comp0(shuffle1.fprim(_b, comp0(_c, _d)), _a),
#             comp0s(shuffle1.fprim(_c, _d), _b, _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             comp0s(_b, shuffle1.fprim(_a, _c), _d),
#             comp0s(_b, _c, shuffle1.fprim(_a, _d)),
#             comp0(shuffle20.fprim(_b, _c, _d), _a),
#             comp0s(shuffle1.fprim(_c, _d), _b, _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             comp0s(_b, shuffle1.fprim(_a, _c), _d),
#             comp0s(shuffle1.fprim(_b, _c), _a, _d),
#             comp0(_c, shuffle21.adj.fprim(_b, _a, _d)),
#             comp0s(shuffle1.fprim(_c, _d), _b, _a)),
#         comp1s(
#             comp0s(shuffle1.fprim(_a, _b), _c, _d),
#             comp0s(_b, shuffle1.fprim(_a, _c), _d),
#             comp0s(shuffle1.fprim(_b, _c), _a, _d),
#             shuffle21.adj.fprim(_c, comp0(_b, _a), _d)),
#         shuffle1Dim2.fprim(comp1s(comp0(shuffle1.fprim(_a, _b), _c), comp0(_b, shuffle1.fprim(_a, _c)), comp0(shuffle1.fprim(_b, _c), _a)), ensureEqMol1(_d)),
#         comp1s(
#             shuffle1.fprim(comp0s(_a, _b, _c), _d),
#             comp0s(_d, shuffle1.fprim(_a, _b), _c),
#             comp0(_d, shuffle21.adj.fprim(_b, _a, _c))),
#         comp1s(
#             shuffle1.fprim(comp0s(_a, _b, _c), _d),
#             comp0(_d, shuffle1Dim2.fprim(shuffle1.fprim(_a, _b), ensureEqMol1(_c))))
#     )
# )


# xx = circleAround3(breen4Target2, [
#     shuffle30,
#     shuffle30.adj,
#     shuffle30.getMate(),
#     shuffle30.adj.getMate(),
#     shuffle31,
#     shuffle31.adj,
#     shuffle31.getMate(),
#     shuffle31.adj.getMate(),
#     shuffle32,
#     shuffle32.adj,
#     shuffle32.getMate(),
#     shuffle32.adj.getMate(),
#     breen3,
#     breen3.adj,
#     breen3.getMate(),
#     breen3.adj.getMate(),
#     shuffle1Dim3,
#     shuffle1Dim3.adj,
#     shuffle1Dim3.getMate(),
#     shuffle1Dim3.adj.getMate(),
#     shuffle20Dim3,
#     shuffle20Dim3.adj,
#     shuffle20Dim3.getMate(),
#     shuffle20Dim3.adj.getMate(),
#     shuffle21Dim3,
#     shuffle21Dim3.adj,
#     shuffle21Dim3.getMate(),
#     shuffle21Dim3.adj.getMate(),
#
#     shuffle20.adjId1,
#     shuffle20.adjId2,
#     shuffle21.adjId1,
#     shuffle21.adjId2,
#     shuffle1Dim2.adjId1,
#     shuffle1Dim2.adjId2])

# xx = findPaths3(breen4Source2, breen4Target2, [
#     shuffle30,
#     shuffle30.adj,
#     shuffle30.getMate(),
#     shuffle30.adj.getMate(),
#     shuffle31,
#     shuffle31.adj,
#     shuffle31.getMate(),
#     shuffle31.adj.getMate(),
#     shuffle32,
#     shuffle32.adj,
#     shuffle32.getMate(),
#     shuffle32.adj.getMate(),
#     breen3,
#     breen3.adj,
#     breen3.getMate(),
#     breen3.adj.getMate(),
#     shuffle1Dim3,
#     shuffle1Dim3.adj,
#     shuffle1Dim3.getMate(),
#     shuffle1Dim3.adj.getMate(),
#     shuffle20Dim3,
#     shuffle20Dim3.adj,
#     shuffle20Dim3.getMate(),
#     shuffle20Dim3.adj.getMate(),
#     shuffle21Dim3,
#     shuffle21Dim3.adj,
#     shuffle21Dim3.getMate(),
#     shuffle21Dim3.adj.getMate(),
#
#     shuffle20.adjId1,
#     shuffle20.adjId2,
#     shuffle21.adjId1,
#     shuffle21.adjId2,
#     shuffle1Dim2.adjId1,
#     shuffle1Dim2.adjId2])
# yy = set(map(primCount, xx))
# print(len(yy))

# def ff(ii):
#     return findPaths3([0], breen4Paths2[ii],
#         [shuffle30, shuffle31, shuffle32.adj, breen3, shuffle1Dim3, shuffle1Dim3,
#         shuffle20Dim3, shuffle20Dim3.adj, shuffle21Dim3, shuffle21Dim3.adj,
#         shuffle20.adjId1, shuffle20.adjId2, shuffle21.adjId1, shuffle21.adjId2, shuffle1Dim2.adjId2, shuffle1Dim2.adjId2,
#         shuffle1Dim3.getMate(), breen3.getMate()])
# breen4Paths3 = ff(1)
# ls = set(map(lambda x: frozenset(primSetCompCell3(x)), breen4Paths3))
# print(len(breen4Paths3))
# print(len(ls))
#print(list(ls)[0])
#print(list(ls)[1])


#exploreFromToAdjCells3(breen4Paths2, breen4Paths2, [shuffle30, shuffle31, shuffle32.adj, breen3.adj, shuffle1Dim3, shuffle20Dim3, shuffle21Dim3], [])
#searchForPathPairs3AdjFam(small, large, [shuffle30, shuffle31, shuffle32.adj, breen3.adj, shuffle1Dim3], [])
#meetInMiddleSearchForPathPairs3(breen4Paths2, [shuffle30, shuffle31, shuffle32, breen3, shuffle1Dim3])
# , , shuffle32, breen3, shuffle1Dim3, shuffle1M, shuffle1Pi

# [(a @ shuffle1(b, c) @ d) . ((a @ c) @ shuffle1(b, d)) . (a @ shuffle1(c, d) @ b) . (shuffle1(a, (d @ c @ b)))] -- [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, c) @ d) . (shuffle1(b, c) @ (a @ d)) . (shuffle1((c @ b @ a), d))]
# [(a @ shuffle1(b, c) @ d) . ((a @ c) @ shuffle1(b, d)) . (a @ shuffle1(c, d) @ b) . (shuffle1(a, (d @ c @ b)))] -- [(a @ shuffle1(b, c) @ d) . (shuffle1(a, c) @ (b @ d)) . (c @ shuffle1(a, b) @ d) . (shuffle1((c @ b @ a), d))]

#meetInMiddleSearchForPathPairs2(breen4Paths1, [shuffle20,  shuffle21.adj, shuffle1Dim2])
