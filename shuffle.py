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
shuffle1 = PrimitiveFamily("shuffle1", 1, [0, 0], shuffle1Source, shuffle1Target)
shuffle1Dim2 = tritrans2CellPrimFamily(shuffle1)
shuffle1Dim3 = tritrans3CellPrimFamily(shuffle1, shuffle1Dim2)
shuffle1Pi = tritransPiPrimFamily(shuffle1, shuffle1Dim2)
shuffle1M = tritransMPrimFamily(shuffle1, shuffle1Dim2)

# We want paths involving only shuffled vertices, so we need to use a custom endFunc
#shuffle20Vertices = list(map(comp0s, shuffles([_a], [_b, _c])))
#shuffle20Paths1 = findPaths(shuffle20Vertices[0], targetEndCond(shuffle20Vertices[-1]), lambda x: x["current"] not in shuffle20Vertices, [shuffle1], cellsAway1, collate1)
# [(shuffle1(a, (b @ c)))] -- [(shuffle1(a, b) @ c) . (b @ shuffle1(a, c))]
shuffle20Source = minimalASTFromEqMol1(ensureEqMol1(
    shuffle1.buildPrim([_a, comp0(_b, _c)])), [_a, _b, _c])
shuffle20Target = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(shuffle1.buildPrim([_a, _b]), _c),
        comp0(_b, shuffle1.buildPrim([_a, _c])))), [_a, _b, _c])
shuffle20 = PrimitiveFamily("shuffle20", 2, [0, 0, 0], shuffle20Source, shuffle20Target)

#shuffle21Vertices = list(map(comp0s, shuffles([_a, _b], [_c])))
#shuffle21Paths1 = findPaths(shuffle21Vertices[0], targetEndCond(shuffle21Vertices[-1]), lambda x: x["current"] not in shuffle21Vertices, [shuffle1], cellsAway1, collate1)
# [(shuffle1((a @ b), c))] -- [(a @ shuffle1(b, c)) . (shuffle1(a, c) @ b)]
shuffle21Source = minimalASTFromEqMol1(ensureEqMol1(
    shuffle1.buildPrim([comp0(_a, _b), _c])), [_a, _b, _c])
shuffle21Target = minimalASTFromEqMol1(ensureEqMol1(
    comp1(
        comp0(_a, shuffle1.buildPrim([_b, _c])),
        comp0(shuffle1.buildPrim([_a, _c]), _b))), [_a, _b, _c])
shuffle21 = PrimitiveFamily("shuffle21", 2, [0, 0, 0], shuffle21Source, shuffle21Target)

# NOTE: For Dim1 in all the below, We picked the codomain candidate that had no tensors inside of the individual shuffles, and the domain candidate that was a single prim.
# ABCD -> BCDA (n = 1, k = 3, insert A into BCD)
#shuffle30Vertices = list(map(comp0s, shuffles([_a], [_b, _c, _d])))
#shuffle30Paths1 = findPaths(shuffle30Vertices[0], targetEndCond(shuffle30Vertices[-1]), lambda x: x["current"] not in shuffle30Vertices, [shuffle1], cellsAway1, collate1)
#searchForPathPairs2(shuffle30Paths1, [shuffle20, shuffle21])
# [(shuffle1(a, (b @ c @ d)))] -- [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle1(a, c) @ d) . ((b @ c) @ shuffle1(a, d))]
#shuffle30Source1 = ensureEqMol1(shuffle1.buildPrim([_a, comp0s([_b, _c, _d])]))
#shuffle30Target1 = ensureEqMol1(comp1s([
#     comp0s([shuffle1.buildPrim([_a, _b]), _c, _d]),
#     comp0s([_b, shuffle1.buildPrim([_a, _c]), _d]),
#     comp0s([_b, _c, shuffle1.buildPrim([_a, _d])])
# ]))
# shuffle30Paths2 = findPaths2(shuffle30Source1, shuffle30Target1, [shuffle20])
# [[shuffle20(a, (b @ c), d)] & [(shuffle20(a, b, c) @ d) . ((b @ c) @ shuffle1(a, d))]]  --  [[shuffle20(a, b, (c @ d))] & [(shuffle1(a, b) @ (c @ d)) . (b @ shuffle20(a, c, d))]]
shuffle30Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        shuffle20.buildPrim([_a, comp0(_b, _c), _d]),
        comp1(
            comp0(shuffle20.buildPrim([_a, _b, _c]), _d),
            comp0s([_b, _c, shuffle1.buildPrim([_a, _d])])
        )
    )), [_a, _b, _c, _d])
shuffle30Target = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2(
        shuffle20.buildPrim([_a, _b, comp0(_c, _d)]),
        comp1(
            comp0s([shuffle1.buildPrim([_a, _b]), _c, _d]),
            comp0(_b, shuffle20.buildPrim([_a, _c, _d]))
        )
    )), [_a, _b, _c, _d])
shuffle30 = PrimitiveFamily("shuffle30", 3, [0, 0, 0, 0], shuffle30Source, shuffle30Target)

# ABCD -> CDAB (n = 2, k = 2, insert AB into CD)
#shuffle31Vertices = list(map(comp0s, shuffles([_a, _b], [_c, _d])))
#shuffle31Paths1 = findPaths(shuffle31Vertices[0], targetEndCond(shuffle31Vertices[-1]), lambda x: x["current"] not in shuffle31Vertices, [shuffle1], cellsAway1, collate1)
#searchForPathPairs2(shuffle31Paths1, [shuffle20, shuffle21])
#[(shuffle1((a @ b), (c @ d)))] -- [(a @ shuffle1(b, c) @ d) . (shuffle1(a, c) @ (b @ d)) . ((c @ a) @ shuffle1(b, d)) . (c @ shuffle1(a, d) @ b)]
shuffle31Source1  = ensureEqMol1(shuffle1.buildPrim([comp0(_a, _b), comp0(_c, _d)]))
shuffle31Target1 = ensureEqMol1(comp1s([
    comp0s([_a, shuffle1.buildPrim([_b, _c]), _d]),
    comp0s([shuffle1.buildPrim([_a,_c]), _b, _d]),
    comp0s([_c, _a, shuffle1.buildPrim([_b, _d])]),
    comp0s([_c, shuffle1.buildPrim([_a, _d]), _b])
]))
shuffle31Paths2 = findPaths2(shuffle31Source1, shuffle31Target1, [shuffle20, shuffle21])
# [[shuffle20((a @ b), c, d)] & [(shuffle21(a, b, c) @ d) . (c @ shuffle1((a @ b), d))] & [(a @ shuffle1(b, c) @ d) . (shuffle1(a, c) @ (b @ d)) . (c @ shuffle21(a, b, d))]]
# --[[shuffle21(a, b, (c @ d))] & [(a @ shuffle20(b, c, d)) . (shuffle1(a, (c @ d)) @ b)] & [(a @ shuffle1(b, c) @ d) . ((a @ c) @ shuffle1(b, d)) . (shuffle20(a, c, d) @ b)]]
shuffle31Source = minimalASTFromEqAEMol2(ensureEqAEMol2(
    comp2s([
    shuffle20.buildPrim([comp0(_a, _b), _c, _d]),
    comp1(
        comp0(shuffle21.buildPrim([_a, _b, _c]), _d),
        comp0(_c, shuffle1.buildPrim([comp0(_a, _b), _d]))
    ),
    comp1s([
        comp0s([_a, shuffle1.buildPrim([_b, _c]), _d]),
        ,
    ])
    ])), [_a, _b, _c, _d])


# shuffle31Dom1 = molecule1 [r|\shuffle1_{(a \t b),(c \t d)}|]
# shuffle31Cod1 = molecule1 [r|(1_{a} \t \shuffle1_{b,c} \t 1_{d}) . (\shuffle1_{a,c} \t 1_{(b \t d)}) . (1_{(c \t a)} \t \shuffle1_{b,d}) . (1_{c} \t \shuffle1_{a,d} \t 1_{b})|]
# shuffle31Dom2 = molecule2 [r|\shuffle20_{(a \t b),c,d} 2. (\shuffle21_{a,b,c} \t 1_{1_{d}}) * (1_{c} \t \shuffle1_{(a \t b),d}) 2. (1_{a} \t \shuffle1_{b,c} \t 1_{d}) . (\shuffle1_{a,c} \t 1_{(b \t d)}) * (1_{1_{c}} \t \shuffle21_{a,b,d})|]
# shuffle31Cod2 = molecule2 [r|\shuffle21_{a,b,(c \t d)} 2. (1_{1_{a}} \t \shuffle20_{b,c,d}) * (\shuffle1_{a,(c \t d)} \t 1_{b}) 2. (1_{a} \t \shuffle1_{b,c} \t 1_{d}) . (1_{(a \t c)} \t \shuffle1_{b,d}) * (\shuffle20_{a,c,d} \t 1_{1_{b}})|]
