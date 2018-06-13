from fill import *
from ontology import *
from primfamily import *
from transfor import *

# Describe tensors with a function on the left and right as primitive families.

_a = ConstPrim0("a")
_b = ConstPrim0("b")
_ap = ConstPrim0("a'")
_bp = ConstPrim0("b'")

_f = ConstPrim1("f", _a, _b)
_g = ConstPrim1("g", _ap, _bp)
_fa = ConstPrim1("fa", _a, _ap)
_fb = ConstPrim1("fb", _b, _bp)
_ga = ConstPrim1("ga", _a, _ap)
_gb = ConstPrim1("gb", _b, _bp)

_ala = ConstPrim2("alpha_a", _fa, _ga)
_alb = ConstPrim2("alpha_b", _fb, _gb)

# TODO: this might be wrong; might need to replace Comp0Node with 
# a PrimitiveFamily node for a new primitive family which represents tensor
# leftTensorId1Source = Comp0Node(
#     SourceNode(VarNode(0)),
#     VarNode(1)
# )
# leftTensorId1Target = Comp0Node(
#     TargetNode(VarNode(0)),
#     VarNode(1)
# )
# leftTensorId1 = PrimitiveFamily("leftTensorId1", 1, 0, [1, 0], leftTensorId1Source, leftTensorId1Target)

# # same for idTensorRight1...
# idTensorRight1Source = Comp0Node(
#     VarNode(0),
#     SourceNode(VarNode(1))
# )
# idTensorRight1Target = Comp0Node(
#     VarNode(0),
#     TargetNode(VarNode(1))
# )
# idTensorRight1 = PrimitiveFamily("idTensorRight1", 1, 0, [0, 1], idTensorRight1Source, idTensorRight1Target)

# # now sigma
# sigma2Source = Comp1Node(
#     PrimitiveFamilyNode(leftTensorId1,[
#         VarNode(0), 
#         SourceNode(VarNode(1))
#     ]),
#     PrimitiveFamilyNode(idTensorRight1,[
#         TargetNode(VarNode(0)), 
#         VarNode(1)
#     ])  
# )

# sigma2Target = Comp1Node(
#     PrimitiveFamilyNode(idTensorRight1,[
#         SourceNode(VarNode(0)), 
#         VarNode(1)
#     ]),
#     PrimitiveFamilyNode(leftTensorId1,[
#         VarNode(0), 
#         TargetNode(VarNode(1))
#     ])
# )

# sigma = PrimitiveFamily("sigma", 2, 0, [1,1], sigma2Source, sigma2Target)

# axioms for sigma: encode as 3-cells
# source, target, primitivefamily







# ALTERNATIVE APPROACH: Using the pseudofunctor B \times B -> B

sigma2Source = Comp0Node(
    VarNode(0),
    VarNode(1)
)

sigma2Target = sigma2Source

sigma = PrimitiveFamily("sigma", 2, 0, [1,1], sigma2Source, sigma2Target)

# In order from K-theory for 2-categories.
# a.alphab sigma alphaa.b
# sigmaAx1Source = comp2s(
#     [
#         comp1(comp0(_a,_alb), comp0(_fa,_bp)),
#         sigma.fprim(_fa,_gb),
#         comp1(comp0(_ala,_b), comp0(_ap,_gb))
#     ]
# )

nala = VarNode(0)
nalb = VarNode(1)
nfa = SourceNode(VarNode(0))
nfb = SourceNode(VarNode(1))
nga = TargetNode(VarNode(0))
ngb = TargetNode(VarNode(1))
na = SourceNode(nfa)
nb = SourceNode(nfb)
nap = TargetNode(nfa)
nbp = TargetNode(nfb)
#_ala is input 0, _alb is input 1 (to the axiom)
sigmaAx1Source = Comp2Node(
    Comp2Node(
        Comp1Node(
            Comp0Node(
                na, nalb
            ),
            Comp0Node(
                nfa, nbp
            )
        ),
        PrimitiveFamilyNode(
            sigma,
            [nfa, ngb]
        )
    ),
    Comp1Node(
        Comp0Node(
            nala, nb
        ),
        Comp0Node(
            nap, ngb
        )
    )
)

sigmaAx1Target = Comp2Node(
    Comp2Node(
        Comp1Node(
            Comp0Node(
                na, nfb
            ),
            Comp0Node(
                nala, nbp
            )
        ),
        PrimitiveFamilyNode(
            sigma,
            [nga, nfb]
        )
    ),
    Comp1Node(
        Comp0Node(
            nga, nb
        ),
        Comp0Node(
            nap, nalb
        )
    )
)

ala = _ala
alb = _alb