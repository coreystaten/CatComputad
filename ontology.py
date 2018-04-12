from enum import Enum


class IdHashed(object):
    def __hash__(self):
        return id(self)

    def __eq__(self,x):
        return id(self) == id(x)


class Prim0(IdHashed):
    pass

class ConstPrim0(Prim0):
    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return self.name

class Mol0(IdHashed):
    def __init__(self, prim0s):
        self.prim0s = prim0s
        self.hash = sum(map(hash,self.prim0s))

    def __str__(self):
        if len(self.prim0s) == 0:
            return "EMPTY-DONT-DISPLAY"
        elif len(self.prim0s) == 1:
            return str(self.prim0s[0])
        else:
            return "(" + " @ ".join(map(str,self.prim0s)) + ")"

    def __repr__(self):
        return "Mol0(" + ", ".join(map(repr, self.prim0s)) + ")"

    def __len__(self):
        return len(self.prim0s)

# NOTE: For equivalence results to work, 1-primitives must not have source or
# target which is not fMol0(()).
class Prim1(IdHashed):
    pass

class ConstPrim1(Prim1):
    def __init__(self,name,source,target):
        self.name = name
        self.source = source
        self.target = target

    def __repr__(self):
        return self.name

class CollapsePrim1(Prim1):
    def __init__(self,prim2):
        self.prim2 = prim2
        self.source = prim2.source.source
        self.target = prim2.target.target

    def __repr__(self):
        return "collapse(" + str(self.prim2) + ")"

class Atom1(IdHashed):
    def __init__(self,a0,p1,b0):
        if not(type(a0) == Mol0 and type(b0) == Mol0 and isinstance(p1, Prim1)):
            raise "ERROR"

        self.a0 = a0
        self.b0 = b0
        self.p1 = p1
        self.source = comp0s(self.a0, self.p1.source, self.b0)
        self.target = comp0s(self.a0, self.p1.target, self.b0)

    def __str__(self):
        part0s = []

        if len(self.a0):
            part0s.append(str(self.a0))
        part0s.append(str(self.p1))
        if len(self.b0):
            part0s.append(str(self.b0))

        if len(part0s) > 1:
            return " @ ".join(part0s)
        else:
            return part0s[0]


    def __repr__(self):
        return "Atom1(" + repr(self.a0) + ", " + repr(self.p1) + ", " + repr(self.b0) + ")"

class Mol1(IdHashed):
    pass

class IdMol1(Mol1):
    def __init__(self, mol0):
        self.mol0 = mol0
        self.atom1s = ()
        self.source = self.mol0
        self.target = self.mol0

    def __str__(self):
        return "1_{" + str(self.mol0) + "}"

    def __repr__(self):
        return "IdMol1(" + repr(self.mol0) + ")"

    def __len__(self):
        return 0


class NonIdMol1(Mol1):
    def __init__(self, atom1s):
        self.atom1s = atom1s
        self.source = self.atom1s[0].source
        self.target = self.atom1s[-1].target

    def __str__(self):
        return " . ".join(map(lambda x: "(" + str(x) + ")", self.atom1s))

    def __repr__(self):
        return "NonIdMol1(" + repr(self.atom1s) + ")"

    def __len__(self):
        return len(self.atom1s)

class EqMol1(IdHashed):
    def __init__(self, mol1s):
        for inst in mol1s:
            if not(isinstance(inst, Mol1)):
                raise Exception("ERROR")

        self.mol1s = mol1s
        self.length = len(next(iter(self.mol1s)))
        self.source = next(iter(self.mol1s)).source
        self.target = next(iter(self.mol1s)).target

    def __str__(self):
        return "[" + str(next(iter(self.mol1s))) + "]"

    def __repr__(self):
        return "EqMol1(" + repr(self.mol1s) + ")"

    def __len__(self):
        return self.length

class Prim2(IdHashed):
    def __init__(self):
        self.collapse = CollapsePrim1(self)


class ConstPrim2(Prim2):
    def __init__(self,name,source,target):
        self.name = name
        self.source = source
        self.target = target
        super(ConstPrim2, self).__init__()

    def __str__(self):
        return "p2(" + self.name + ")"

    def __repr__(self):
        return self.name

class Atom2(IdHashed):
    def __init__(self,l1,a0,p2,b0,r1):
        if not(isinstance(l1, Mol1) and isinstance(a0, Mol0) and isinstance(p2, Prim2) and isinstance(b0, Mol0) and isinstance(r1, Mol1)):
            raise Exception("Invalid Atom2 construction.")

        self.l1 = l1
        self.a0 = a0
        self.p2 = p2
        self.b0 = b0
        self.r1 = r1
        self.source = comp1s(self.l1, comp0s(self.a0, self.p2.source, self.b0), self.r1)
        if not(isinstance(self.source, Mol1)):
            raise Exception("Invalid source")
        self.target = comp1s(self.l1, comp0s(self.a0, self.p2.target, self.b0), self.r1)
        self.collapse = comp1s(self.l1, comp0s(self.a0, self.p2.collapse, self.b0), self.r1)

    def __str__(self):
        part1s = []
        part0s = []
        if len(self.l1):
            part1s.append(str(self.l1))
        if len(self.a0):
            part0s.append(str(self.a0))
        part0s.append(str(self.p2))
        if len(self.b0):
            part0s.append(str(self.b0))
        if len(part0s) > 1:
            part1s.append("(" + " @ ".join(part0s) + ")")
        else:
            part1s.append(part0s[0])
        if len(self.r1):
            part1s.append(str(self.r1))
        if len(part1s) > 1:
            return " . ".join(part1s)
        else:
            return part1s[0]


    def __repr__(self):
        return "Atom2(" + repr(self.l1) + ", " + repr(self.a0) + ", " + repr(self.p2) + ", " + repr(self.b0) + ", " + repr(self.r1) + ")"

class EqAtom2(IdHashed):
    def __init__(self,atom2s):
        for inst in atom2s:
            if not(isinstance(inst, Atom2)):
                raise Exception("ERROR")

        self.atom2s = atom2s
        self.righthand = min(self.atom2s, key=lambda x: len(x.l1.atom1s))
        self.lefthand = min(self.atom2s, key=lambda x: len(x.r1.atom1s))
        self.source = fEqMol1(next(iter(self.atom2s)).source)
        self.target = fEqMol1(next(iter(self.atom2s)).target)

    def __str__(self):
        return "[" + str(self.righthand) + "]"

    def __repr__(self):
        return "EqAtom2(" + repr(self.atom2s) + ")"

class AEMol2(IdHashed):
    pass

class AEIdMol2(AEMol2):
    def __init__(self, eqMol1):
        if type(eqMol1) != EqMol1:
            raise "ERROR"

        self.eqMol1 = eqMol1
        self.eqAtom2s = ()
        self.source = self.eqMol1
        self.target = self.eqMol1

    def __str__(self):
        return "1_{" + str(self.eqMol1) + "}"

    def __repr__(self):
        return "AEIdMol2(" + repr(self.eqMol1) + ")"

    def __len__(self):
        return 0

class AENonIdMol2(AEMol2):
    def __init__(self, eqAtom2s):
        self.eqAtom2s = eqAtom2s
        self.source = self.eqAtom2s[0].source
        self.target = self.eqAtom2s[-1].target

    def __str__(self):
        if len(self.eqAtom2s) > 1:
            return " & ".join(map(str, self.eqAtom2s))
        else:
            return str(self.eqAtom2s[0])

    def __repr__(self):
        return "AEMol2(" + repr(list(self.eqAtom2s)) + ")"

    def __len__(self):
        return len(self.eqAtom2s)

class EqAEMol2(IdHashed):
    def __init__(self, aeMol2s):
        self.aeMol2s = aeMol2s
        self.length = len(next(iter(self.aeMol2s)))
        self.source = next(iter(aeMol2s)).source
        self.target = next(iter(aeMol2s)).target

    def __str__(self):
        return "[" + str(next(iter(self.aeMol2s))) + "]"

    def __repr__(self):
        return "EqAEMol2(" + repr(list(self.aeMol2s)) + ")"

    def __len__(self):
        return self.length

class Prim3(IdHashed):
    pass

class ConstPrim3(Prim3):
    def __init__(self,name,source,target):
        self.name = name
        self.source = source
        self.target = target

    def __str__(self):
        return "p3(" + self.name + ")"

    def __repr__(self):
        return self.name

# Components of StepCell3 are up to equivalence, since they're just an output format.
# Composition order is the same as for a 3-atom.
class StepCell3(IdHashed):
    def __init__(top, l, a, p, b, r, bot):
        self.top = top
        self.l = l
        self.a = a
        self.p = p
        self.b = b
        self.r = r
        self.bot = bot
        self.source = comp2s(top, comp1s(l, comp0s(a, p3.source, b), r), bot)
        self.target = comp2s(top, comp1s(l, comp0s(a, p3.target, b), r), bot)

    def __str__(self):
        _top = asLowestDim(top)
        _l = asLowestDim(l)
        _a = asLowestDim(a)
        _b = asLowestDim(b)
        _r = asLowestDim(r)
        _bot = asLowestDim(bot)

        part2s = []
        part1s = []
        part0s = []

        if not(isinstance(_top, EqMol1)):
            part2s.append(str(_top))
        if not(isinstance(_l, Mol0)):
            part1s.append(str(_l))
        if not(isinstance(_a, Mol0) and len(a) == 0):
            part0s.append(str(_a))
        part0s.append(str(self.p))
        if not(isinstance(_b, Mol0) and len(b) == 0):
            part0s.append(str(_b))
        if len(part0s) > 1:
            part1s.append("(" + " @ ".join(part0s) + ")")
        else:
            part1s.append(part0s[0])
        if not(isinstance(_r, Mol0)):
            part1s.append(str(_r))
        if len(part1s) > 1:
            part2s.append("(" + " . ".join(part1s) + ")")
        else:
            part2s.append(part1s[0])
        if not(isinstance(_bot, EqMol1)):
            part2s.append(str(_bot))
        if len(part2s) > 1:
            return " & ".join(part2s)
        else:
            return part2s[0]


    def __repr__(self):
        return "StepCell3(" + repr(self.top) + ", " + repr(self.l) + ", " + repr(self.a) + ", " + repr(self.p) + ", " + repr(self.b) + ", " + repr(self.r) + ", " + repr(self.bot) + ")"

class CompCell3(IdHashed):
    def __init__(self, stepCell3s):
        self.stepCell3s = stepCell3s
        self.source = self.stepCell3s[0].source
        self.target = self.stepCell3s[-1].target

    def __len__(self):
        return len(self.stepCell3s)

    def __repr__(self):
        return "CompCell3(" + ", ".join(map(repr, self.stepCell3s)) + ")"

class Functor(IdHashed):
    def __init__(self, name):
        self.name = name

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Functor(" + self.name + ")"

# FunctorPrims represent cells inside of a functor application.
class FunctorPrim0(Prim0):
    def __init__(self, functor, mol0):
        self.functor = functor
        self.mol0 = mol0

    def __str__(self):
        return str(self.functor) + "(" + str(self.mol0) + ")"

    def __repr__(self):
        return "FunctorPrim0(" + repr(self.functor) + ", " + repr(self.mol0) + ")"

class FunctorPrim1(Prim1):
    def __init__(self, functor, eqMol1):
        self.functor = functor
        self.eqMol1 = eqMol1
        self.source = fFunctorPrim0(functor, eqMol1.source)
        self.target = fFunctorPrim0(functor, eqMol1.target)

    def __str__(self):
        return str(self.functor) + "(" + str(self.eqMol1) + ")"

    def __repr__(self):
        return "FunctorPrim1(" + repr(self.functor) + ", " + repr(eself.qMol1) + ")"

class FunctorPrim2(Prim2):
    def __init__(self, functor, eqAEMol2):
        self.functor = functor
        self.eqAEMol2 = eqAEMol2
        self.source = fFunctorPrim1(functor, eqAEMol2.source)
        self.target = fFunctorPrim1(functor, eqAEMol2.target)

    def __str__(self):
        return str(self.functor) + "(" + str(self.eqAEMol2) + ")"

    def __repr__(self):
        return "FunctorPrim2(" + repr(self.functor) + ", " + repr(self.eqAEMol2) + ")"

class FunctorPrim3(Prim3):
    def __init__(self, functor, cell3):
        self.functor = functor
        self.cell3 = cell3
        self.source = fFunctorPrim2(functor, cell3.source)
        self.target = fFunctorPrim2(functor, cell3.target)

    def __str__(self):
        return str(self.functor) + "(" + str(self.cell3) + ")"

    def __repr__(self):
        return "FunctorPrim3(" + repr(self.functor) + ", " + repr(self.cell3) + ")"


class TransposeType(Enum):
    LEFT = 0
    RIGHT = 1
    NONE = 2

def take(l, n):
    if n >= len(l):
        return l
    elif n <= 0:
        return ()
    else:
        return l[:n]

def drop(l, n):
    if n >= len(l):
        return ()
    elif n <= 0:
        return l
    else:
        return l[n:]

def takeMol0(mol0, n):
    if n <= 0:
        return fMol0(())
    return fMol0(take(mol1.prim0s, n))

def dropMol0(mol0, n):
    if n >= len(mol0):
        return fMol0(())
    return fMol0(drop(mol1.prim0s, n))

def takeMol1(mol1, n):
    if n <= 0:
        return fIdMol1(mol1.source)
    return fNonIdMol1(take(mol1.atom1s, n))

def dropMol1(mol1, n):
    if n >= len(mol1):
        return fIdMol1(mol1.target)
    return fNonIdMol1(drop(mol1.atom1s, n))

def takeAEMol2(aeMol2, n):
    if n <= 0:
        return fAEIdMol2(aeMol2.source)
    return fAENonIdMol2(take(aeMol2.eqAtom2s, n))

def dropAEMol2(aeMol2, n):
    if n >= len(aeMol2):
        return fAEIdMol2(aeMol2.target)
    return fAENonIdMol2(drop(aeMol2.eqAtom2s, n))

def prim0ToMol0(p):
    return fMol0((p,))

def prim1ToAtom1(p):
    return fAtom1(fMol0(()), p, fMol0(()))

def prim1ToMol1(p):
    return fNonIdMol1((prim1ToAtom1(p),))

def atom1ToMol1(a):
    return fNonIdMol1((a,))

def prim2ToAtom2(p):
    return fAtom2(fIdMol1(p.source.source), fMol0(()), p, fMol0(()), fIdMol1(p.target.target))

def prim2ToEqAtom2(p):
    return fEqAtom2(prim2ToAtom2(p))

def prim2ToAEMol2(p):
    return fAENonIdMol2((prim2ToEqAtom2(p),))

# asPrimN functions take an equivalence cell (Mol0, EqMol1, or EqAEMol2) and attempt to cast them as a single primitive.
# If they are not just a single primitive, they return None.
def asPrim0(mol0):
    if len(mol0) == 1:
        return mol0.prim0s[0]
    return None

def asPrim1(eqMol1):
    inst = next(iter(eqMol1.mol1s))
    if len(inst) == 1:
        atom = inst.atom1s[0]
        if len(atom.a0) == 0 and len(atom.b0) == 0:
            return atom.p1
    return None

def asPrim2(eqAEMol2):
    inst = next(iter(eqAEMol2.aeMol2s))
    if len(inst) == 1:
        eqAtom2 = inst.eqAtom2s[0]
        atInst = next(iter(eqAtom2.atom2s))
        if len(atInst.a0) == 0 and len(atInst.b0) == 0 and len(atInst.l1) == 0 and len(atInst.r1) == 0:
            return atom.p2
    return None

# Takes any cell and, if it's an identity, returns its lowest dimensional representative.  Otherwise returns the cell.
def asLowestDim(cell):
    result = cell
    if isinstance(cell, EqAEMol2):
        cell = next(iter(cell.aeMol2s))
    if isinstance(cell, AEIdMol2):
        cell = cell.eqMol1
        result = cell
    if isinstance(cell, EqMol1):
        cell = next(iter(cell.mol1s))
    if isinstance(cell, IdMol1):
        cell = cell.mol0
        result = cell
    return result

def isIdMol1(p):
    if isinstance(p, IdMol1) or (isinstance(p, EqMol1) and isinstance(next(iter(p.mol1s)), IdMol1)):
        return True

def source0(x):
    if type(x) in [Atom1, IdMol1, NonIdMol1, EqMol1] or isinstance(x, Prim1):
        return x.source
    if type(x) in [Atom2, EqAtom2, AENonIdMol2, AEIdMol2, EqAEMol2] or isinstance(x, Prim2):
        return x.source.source
    raise Exception("Source0 of unknown type %s" % (type(x),))

def target0(x):
    if type(x) in [Atom1, IdMol1, NonIdMol1, EqMol1] or isinstance(x, Prim1):
        return x.target
    if type(x) in [Atom2, EqAtom2, AENonIdMol2, AEIdMol2, EqAEMol2] or isinstance(x, Prim2):
        return x.target.target
    raise Exception("Target0 of unknown type %s" % (type(x),))

comp0Table = {
  (Prim0, Prim0): lambda x, y: comp0(prim0ToMol0(x), prim0ToMol0(y)),
  (Mol0, Mol0): lambda x, y: fMol0(x.prim0s + y.prim0s),

  (Mol0, Prim1): lambda x, y: comp0(x, prim1ToAtom1(y)),
  (Prim1, Mol0): lambda x, y: comp0(prim1ToAtom1(x), y),
  (Mol0, Atom1): lambda x, y: fAtom1(comp0(x, y.a0), y.p1, y.b0),
  (Atom1, Mol0): lambda x, y: fAtom1(x.a0, x.p1, comp0(x.b0, y)),
  (Mol0, NonIdMol1): lambda x, y: fNonIdMol1(tuple(map(lambda z: comp0(x, z), y.atom1s))),
  (NonIdMol1, Mol0): lambda x, y: fNonIdMol1(tuple(map(lambda z: comp0(z, y), x.atom1s))),
  (Mol0, IdMol1): lambda x, y: fIdMol1(comp0(x, y.mol0)),
  (IdMol1, Mol0): lambda x, y: fIdMol1(comp0(x.mol0, y)),
  (IdMol1, IdMol1): lambda x, y: fIdMol1(comp0(x.mol0, y.mol0)),
  (IdMol1, NonIdMol1): lambda x, y: fNonIdMol1(tuple(map(lambda z: comp0(x.mol0, z), y.atom1s))),
  (NonIdMol1, IdMol1): lambda x, y: fNonIdMol1(tuple(map(lambda z: comp0(z, y.mol0), x.atom1s))),
  (NonIdMol1, NonIdMol1): lambda x, y: fNonIdMol1(tuple(map(lambda z: comp0(z, y.source), x.atom1s)) + tuple(map(lambda z: comp0(x.target, z), y.atom1s))),

  (Mol0, Prim2): lambda x, y: comp0(x, prim2ToAtom2(y)),
  (Prim2, Mol0): lambda x, y: comp0(prim2ToAtom2(x), y),
  (Mol0, Atom2): lambda x, y: fAtom2(comp0(x, y.l1), comp0(x, y.a0), y.p2, y.b0, comp0(x, y.r1)),
  (Atom2, Mol0): lambda x, y: fAtom2(comp0(x.l1, y), x.a0, x.p2, comp0(x.b0, y), comp0(x.r1, y)),
  (EqMol1, EqAtom2): lambda x, y: fEqAtom2(comp0(next(iter(x.mol1s)), next(iter(y.atom2s)))),
  (EqAtom2, EqMol1): lambda x, y: fEqAtom2(comp0(next(iter(x.atom2s)), next(iter(y.mol1s)))),
  (NonIdMol1, EqAtom2): lambda x, y: fEqAtom2(comp0(x, next(iter(y.atom2s)))),
  (EqAtom2, NonIdMol1): lambda x, y: fEqAtom2(comp0(next(iter(x.atom2s)), y)),
  (IdMol1, EqAtom2): lambda x, y: fEqAtom2(comp0(x, next(iter(y.atom2s)))),
  (EqAtom2, IdMol1): lambda x, y: fEqAtom2(comp0(next(iter(x.atom2s)), y)),
  (NonIdMol1, Atom2): lambda x, y: fAtom2(comp0(x, y.l1), comp0(x.target, y.a0), y.p2, y.b0, comp0(x.target, y.r1)),
  (Atom2, NonIdMol1): lambda x, y: fAtom2(comp0(x.l1, y.source), x.a0, x.p2, comp0(x.b0, y.source), comp0(x.r1, y)),
  (IdMol1, Atom2): lambda x, y: fAtom2(comp0(x, y.l1), comp0(x.target, y.a0), y.p2, y.b0, comp0(x.target, y.r1)),
  (Atom2, IdMol1): lambda x, y: fAtom2(comp0(x.l1, y.source), x.a0, x.p2, comp0(x.b0, y.source), comp0(x.r1, y)),
  (AENonIdMol2, AENonIdMol2): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp0(z, y.source), x.eqAtom2s)) + tuple(map(lambda z: comp0(x.target, z), y.eqAtom2s)))
}

dimByType = {
    Mol0: 0,
    Atom1: 1,
    IdMol1: 1,
    NonIdMol1: 1,
    EqMol1: 1,
    Atom2: 2,
    EqAtom2: 2,
    AEIdMol2: 2,
    AENonIdMol2: 2,
    EqAEMol2: 2,
    Prim3: 3,
    StepCell3: 3,
    CompCell3: 3
}

def dim(x):
    if isinstance(x, Prim0):
        return 0
    elif isinstance(x, Prim1):
        return 1
    elif isinstance(x, Prim2):
        return 2
    elif isinstance(x, Prim3):
        return 3
    else:
        return dimByType[type(x)]




def comp0(x,y):
    maxDim = max(dim(x), dim(y))

    if maxDim == 3:
        # Handle this case special, because we don't want to escalate things to equivalence classes.
        raise Exception("Unimplemented")

    if maxDim >= 0:
        if isinstance(x, Prim0):
            x = fMol0((x,))
        if isinstance(y, Prim0):
            y = fMol0((y,))

        if isinstance(x, Mol0) and isinstance(y, Mol0):
            return fMol0(x.prim0s + y.prim0s)
    if maxDim >= 1:
        if isinstance(x, Mol0):
            x = fIdMol1(x)
        if isinstance(y, Mol0):
            y = fIdMol1(y)

        if isinstance(x, Prim1):
            x = prim1ToAtom1(x)
        if isinstance(y, Prim1):
            y = prim1ToAtom1(y)

        if isinstance(x, IdMol1) and isinstance(y, IdMol1):
            return fIdMol1(comp0(x.mol0, y.mol0))
        elif isinstance(x, IdMol1) and isinstance(y, Atom1):
            return fAtom1(comp0(x.mol0, y.a0), y.p1, y.b0)
        elif isinstance(x, Atom1) and isinstance(y, IdMol1):
            return fAtom1(x.a0, x.p1, comp0(x.b0, y.mol0))
        elif isinstance(x, IdMol1) and isinstance(y, NonIdMol1):
            return fNonIdMol1(tuple(map(lambda z: comp0(x, z), y.atom1s)))
        elif isinstance(x, NonIdMol1) and isinstance(y, IdMol1):
            return fNonIdMol1(tuple(map(lambda z: comp0(z, y), x.atom1s)))

        if isinstance(x, Atom1):
            x = fNonIdMol1((x,))
        if isinstance(y, Atom1):
            y = fNonIdMol1((y,))

        if isinstance(x, NonIdMol1) and isinstance(y, NonIdMol1):
            return fNonIdMol1(tuple(map(lambda z: comp0(z, y.source), x.atom1s)) + tuple(map(lambda z: comp0(x.target, z), y.atom1s)))

        if isinstance(x, Mol1) and isinstance(y, EqMol1):
            return fEqMol1(comp0(x, next(iter(y.mol1s))))
        elif isinstance(x, EqMol1) and isinstance(y, Mol1):
            return fEqMol1(comp0(next(iter(x.mol1s)), y))
        elif isinstance(x, EqMol1) and isinstance(y, EqMol1):
            return fEqMol1(comp0(next(iter(x.mol1s)), next(iter(y.mol1s))))

    if maxDim >= 2:
        if isinstance(x, Prim2):
            x = prim2ToAtom2(x)
        if isinstance(y, Prim2):
            y = prim2ToAtom2(y)

        if isinstance(x, Atom2) and isinstance(y, Mol1):
            return fAtom2(comp0(x.l1, y.source), x.a0, x.p2, comp0(x.b0, y.source), comp0(x.r1, y))
        elif isinstance(x, Mol1) and isinstance (y, Atom2):
            return fAtom2(comp0(x, y.l1), comp0(x.target, y.a0), y.p2, y.b0, comp0(x.target, y.r1))

        if isinstance(x, Mol1):
            x = fAEIdMol2(fEqMol1(x))
        elif isinstance(x, EqMol1):
            x = fAEIdMol2(x)
        elif isinstance(x, Atom2):
            x = fEqAtom2(x)
        if isinstance(y, Mol1):
            y = fAEIdMol2(fEqMol1(y))
        elif isinstance(y, EqMol1):
            y = fAEIdMol2(y)
        elif isinstance(y, Atom2):
            y = fEqAtom2(y)

        if isinstance(x, AEIdMol2) and isinstance(y, AEIdMol2):
            return fAEIdMol2(comp0(x.eqMol1, y.eqMol1))
        elif isinstance(x, EqAtom2) and isinstance(y, AEIdMol2):
            return fEqAtom2(comp0(next(iter(x.atom2s)), next(iter(y.eqMol1.mol1s))))
        elif isinstance(x, AEIdMol2) and isinstance (y, EqAtom2):
            return fEqAtom2(comp0(next(iter(x.eqMol1.mol1s)), next(iter(y.atom2s))))
        elif isinstance(x, AENonIdMol2) and isinstance(y, AEIdMol2):
            return fAENonIdMol2(tuple(map(lambda z: comp0(z, y), x.eqAtom2s)))
        elif isinstance(x, AEIdMol2) and isinstance(y, AENonIdMol2):
            return fAENonIdMol2(tuple(map(lambda z: comp0(x, z), y.eqAtom2s)))

        if isinstance(x, EqAtom2):
            x = fAENonIdMol2((x,))
        if isinstance(y, EqAtom2):
            y = fAENonIdMol2((y,))

        if isinstance(x, AENonIdMol2) and isinstance(y, AENonIdMol2):
            return fAENonIdMol2(tuple(map(lambda z: comp0(z, y.source), x.eqAtom2s)) + tuple(map(lambda z: comp0(x.target, z), y.eqAtom2s)))

        if isinstance(x, AEMol2) and isinstance(y, EqAEMol2):
            return fEqAEMol2(comp0(x, next(iter(y.aeMol2s))))
        elif isinstance(x, EqAEMol2) and isinstance(y, AEMol2):
            return fEqAEMol2(comp0(next(iter(x.aeMol2s)), y))
        elif isinstance(x, EqAEMol2) and isinstance(y, EqAEMol2):
            return fEqAEMol2(comp0(next(iter(x.aeMol2s)), next(iter(y.aeMol2s))))

    raise Exception("Undetermined comp0 of types %s and %s" % (str(type(x)), str(type(y))))











#def comp0(x,y):
#    if isinstance(x, Prim0):
#        typeX = Prim0
#    elif isinstance(x, Prim1):
#        typeX = Prim1
#    elif isinstance(x, Prim2):
#        typeX = Prim2
#    else:
#        typeX = type(x)
#    if isinstance(y, Prim0):
#        typeY = Prim0
#    elif isinstance(y, Prim1):
#        typeY = Prim1
#    elif isinstance(y, Prim2):
#        typeY = Prim2
#    else:
#        typeY = type(y)
#    typePair = (typeX, typeY)
#
#    if typePair in comp0Table:
#        return comp0Table[typePair](x,y)
#    else:
#        raise Exception("Invalid 0-composition between %s and %s" % (str(type(x)), str(type(y))))

def comp0s(*xs):
    result = xs[0]
    for x in xs[1:]:
        result = comp0(result,x)
    return result

comp1Table = {
    (IdMol1, Atom1): lambda x, y: fNonIdMol1((y,)),
    (Atom1, IdMol1): lambda x, y: fNonIdMol1((x,)),
    (NonIdMol1, Atom1): lambda x, y: fNonIdMol1(x.atom1s + (y,)),
    (Atom1, NonIdMol1): lambda x, y: fNonIdMol1((x,) + y.atom1s),
    (IdMol1, IdMol1): lambda x, y: x,
    (NonIdMol1, IdMol1): lambda x, y: x,
    (IdMol1, NonIdMol1): lambda x, y: y,
    (NonIdMol1, NonIdMol1): lambda x, y: fNonIdMol1(x.atom1s + y.atom1s),
    (EqMol1, EqMol1): lambda x, y: fEqMol1(comp1(next(iter(x.mol1s)), next(iter(y.mol1s)))),

    (IdMol1, Atom2): lambda x, y: y,
    (Atom2, IdMol1): lambda x, y: x,
    (NonIdMol1, Atom2): lambda x, y: fAtom2(comp1(x, y.l1), y.a0, y.p2, y.b0, y.r1),
    (Atom2, NonIdMol1): lambda x, y: fAtom2(x.l1, x.a0, x.p2, x.b0, comp1(x.r1, y)),
    (IdMol1, EqAtom2): lambda x, y: y,
    (EqAtom2, IdMol1): lambda x, y: x,
    (NonIdMol1, EqAtom2): lambda x, y: fEqAtom2(comp1(x, next(iter(y.atom2s)))),
    (EqAtom2, NonIdMol1): lambda x, y: fEqAtom2(comp1(next(iter(x.atom2s)), y)),
    (EqAtom2, EqMol1): lambda x, y: fEqAtom2(comp1(next(iter(x.atom2s)),next(iter(y.mol1s)))),
    (EqMol1, EqAtom2): lambda x, y: fEqAtom2(comp1(next(iter(x.mol1s)),next(iter(y.atom2s)))),

    (AENonIdMol2, EqMol1): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp1(z, y), x.eqAtom2s))),
    (EqMol1, AENonIdMol2): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp1(x, z), y.eqAtom2s))),
    (NonIdMol1, AENonIdMol2): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp1(x, z), y.eqAtom2s))),
    (AENonIdMol2, NonIdMol1): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp1(z, y), x.eqAtom2s))),
    (IdMol1, AENonIdMol2): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp1(x, z), y.eqAtom2s))),
    (AENonIdMol2, IdMol1): lambda x, y: fAENonIdMol2(tuple(map(lambda z: comp1(z, y), x.eqAtom2s))),
    (AENonIdMol2, AENonIdMol2): lambda x, y: comp2(comp1(x, y.source), comp1(x.target, y))
}

def comp1(x,y):
    if target0(x) != source0(y):
        raise Exception("Source target mismatch in 1-composition: %s and %s" % (repr(x), repr(y)))

    typeX = type(x)
    typeY = type(y)
    maxDim = max(dim(x), dim(y))

    if maxDim == 3:
        raise Exception("Unimplemented")

    if maxDim >= 1:
        if isinstance(x, Prim1):
            x = prim1ToMol1(x)
        elif isinstance(x, Atom1):
            x = fNonIdMol1((x,))
        if isinstance(y, Prim1):
            y = prim1ToMol1(y)
        elif isinstance(y, Atom1):
            y = fNonIdMol1((y,))

        if isinstance(x, IdMol1) and isinstance(x, Mol1):
            return y
        elif isinstance(x, Mol1) and isinstance(y, IdMol1):
            return x
        elif isinstance(x, NonIdMol1) and isinstance(y, NonIdMol1):
            return fNonIdMol1(x.atom1s + y.atom1s)

        if isinstance(x, EqMol1) and isinstance(y, EqMol1):
            return fEqMol1(comp1(next(iter(x.mol1s)), next(iter(y.mol1s))))
        elif isinstance(x, Mol1) and isinstance(y, EqMol1):
            return fEqMol1(comp1(x, next(iter(y.mol1s))))
        elif isinstance(x, EqMol1) and isinstance(y, Mol1):
            return fEqMol1(comp1(next(iter(x.mol1s)), y))
    if maxDim == 2:
        if isinstance(x, Mol1) and isinstance(y, Atom2):
            return fAtom2(comp1(x, y.l1), y.a0, y.p2, y.b0, y.r1)
        elif isinstance(x, Atom2) and isinstance(y, Mol1):
            return fAtom2(x.l1, x.a0, x.p2, x.b0, comp1(x.r1, y))

        if isinstance(x, Mol1):
            x = fAEIdMol2(fEqMol1(x))
        elif isinstance(x, EqMol1):
            x = fAEIdMol2(x)
        if isinstance(y, Mol1):
            y = fAEIdMol2(fEqMol1(y))
        elif isinstance(y, EqMol1):
            y = fAEIdMol2(y)

        if isinstance(x, AEIdMol2) and isinstance(y, EqAtom2):
            inst = y.righthand
            return fEqAtom2(fAtom2(comp1(next(iter(x.eqMol1.mol1s)), inst.l1), inst.a0, inst.p2, inst.b0, inst.r1))
        elif isinstance(x, EqAtom2) and isinstance(y, AEIdMol2):
            inst = x.righthand
            return fEqAtom2(fAtom2(inst.l1, inst.a0, inst.p2, inst.b0, comp1(inst.r1, next(iter(y.eqMol1.mol1s)))))

        if isinstance(x, Atom2):
            x = fAENonIdMol2((fEqAtom2(x),))
        elif isinstance(x, EqAtom2):
            x = fAENonIdMol2((x,))
        if isinstance(y, Atom2):
            y = fAENonIdMol2((fEqAtom2(y),))
        elif isinstance(y, EqAtom2):
            y = fAENonIdMol2((y,))

        if isinstance(x, AEIdMol2) and isinstance(y, AEIdMol2):
            return fAEIdMol2(comp1(x.mol1, y.mol1))
        elif isinstance(x, AEIdMol2) and isinstance(y, AENonIdMol2):
            return fAENonIdMol2(tuple(map(lambda z: comp1(x, z), y.eqAtom2s)))
        elif isinstance(x, AENonIdMol2) and isinstance(y, AEIdMol2):
            return fAENonIdMol2(tuple(map(lambda z: comp1(z, y), x.eqAtom2s)))
        elif isinstance(x, AENonIdMol2) and isinstance(y, AENonIdMol2):
            return comp2(comp1(x, y.source), comp1(x.target, y))

        if isinstance(x, AEMol2):
            x = fEqAEMol2(x)
        if isinstance(y, AEMol2):
            y = fEqAEMol2(y)

        if isinstance(x, EqAEMol2) and isinstance(y, EqAEMol2):
            return fEqAEMol2(comp1(next(iter(x.aeMol2s)), next(iter(y.aeMol2s))))

    raise Exception("Undetermined comp1 of types %s and %s" % (str(typeX), str(typeY)))


def comp1s(*xs):
    result = xs[0]
    for x in xs[1:]:
        result = comp1(result,x)
    return result

comp2Table = {
    (AEIdMol2, AEIdMol2): lambda x, y: x,
    (AENonIdMol2, AEIdMol2): lambda x, y: x,
    (AEIdMol2, AENonIdMol2): lambda x, y: y,
    (AENonIdMol2, AENonIdMol2): lambda x, y: fAENonIdMol2(x.eqAtom2s + y.eqAtom2s)
}


def comp2(x,y):
    if ensureEqMol1(x.target) != ensureEqMol1(y.source):
        raise Exception("Source target mismatch in 2-composition: %s and %s" % (repr(x), repr(y)))

    maxDim = max(dim(x), dim(y))

    if maxDim == 2:
        if isinstance(x, Prim2):
            x = prim2ToAEMol2(x)
        elif isinstance(x, Atom2):
            x = fAENonIdMol2((fEqAtom2(x),))
        elif isinstance(x, EqAtom2):
            x = fAENonIdMol2((x,))
        if isinstance(y, Prim2):
            y = prim2ToAEMol2(y)
        elif isinstance(y, Atom2):
            y = fAENonIdMol2((fEqAtom2(y),))
        elif isinstance(y, EqAtom2):
            y = fAENonIdMol2((y,))

        if isinstance(x, AEIdMol2) and isinstance(y, AEMol2):
            return y
        elif isinstance(x, AEMol2) and isinstance(y, AEIdMol2):
            return x
        elif isinstance(x, AENonIdMol2) and isinstance(y, AENonIdMol2):
            return fAENonIdMol2(x.eqAtom2s + y.eqAtom2s)
        elif isinstance(x, EqAEMol2) and isinstance(y, AEMol2):
            return fEqAEMol2(comp2(next(iter(x.aeMol2s)), y))
        elif isinstance(x, AEMol2) and isinstance(y, EqAEMol2):
            return fEqAEMol2(comp2(x, next(iter(y.aeMol2s))))
        elif isinstance(x, EqAEMol2) and isinstance(y, EqAEMol2):
            return fEqAEMol2(comp2(next(iter(x.aeMol2s)), next(iter(y.aeMol2s))))


    raise Exception("Undetermined comp2 of types %s and %s" % (str(type(x)), str(type(y))))


def comp2s(*xs):
    result = xs[0]
    for x in xs[1:]:
        result = comp2(result,x)
    return result

def comp3(x,y):
    if isinstance(x, StepCell3):
        x = CompCell3([x])
    if isinstance(y, StepCell3):
        y = CompCell3([y])

    return CompCell3(x.stepCell3s + y.stepCell3s)


def buildKeyedBy(repo, construct):
    def buildFunc(*args):
        k = tuple(args)

        if k in repo:
            return repo[k]
        else:
            obj = construct(*k)
            repo[k] = obj
            return obj
    return buildFunc

fFunctorPrim0Repo = {}
fFunctorPrim0 = buildKeyedBy(fFunctorPrim0Repo, FunctorPrim0)

fFunctorPrim1Repo = {}
fFunctorPrim1 = buildKeyedBy(fFunctorPrim1Repo, FunctorPrim1)

fFunctorPrim2Repo = {}
fFunctorPrim2 = buildKeyedBy(fFunctorPrim2Repo, FunctorPrim2)

fFunctorPrim3Repo = {}
fFunctorPrim3 = buildKeyedBy(fFunctorPrim3Repo, FunctorPrim3)

mol0Repo = {}
fMol0 = buildKeyedBy(mol0Repo, Mol0)

atom1Repo = {}
fAtom1 = buildKeyedBy(atom1Repo, Atom1)
idMol1Repo = {}
fIdMol1 = buildKeyedBy(idMol1Repo, IdMol1)
nonIdMol1Repo = {}
fNonIdMol1 = buildKeyedBy(nonIdMol1Repo, NonIdMol1)


# If prefix is a prefix of mol0, returns mol0 with the prefix removed.  Otherwise, returns None.
def removePrefixMol0(prefix,mol0):
    p = prefix.prim0s
    m = mol0.prim0s
    if m[:len(p)] == p:
        return fMol0(m[len(p):])
    else:
        return None

# If suffix is a suffix of mol0, returns mol0 with the suffix removed.  Otherwise, returns None.
def removeSuffixMol0(suffix, mol0):
    s = suffix.prim0s
    m = mol0.prim0s
    if m[-len(s):] == s:
        return fMol0(m[:-len(s)])
    else:
        return None

# Returns a tuple of (a2',a1',transposeType), where a2' and a1' may be None.
# By non-degeneracy, they are either left transposable or right, not both.
def transposeAtom1s(a1,a2):
    xLeft1 = removeSuffixMol0(comp0(a2.p1.source, a2.b0), a1.b0)
    xLeft2 = removePrefixMol0(comp0(a1.a0, a1.p1.target), a2.a0)

    if xLeft1 is not None and xLeft1 == xLeft2:
        a2New = fAtom1(comp0s(a1.a0, a1.p1.source, xLeft1), a2.p1, a2.b0)
        a1New = fAtom1(a1.a0, a1.p1, comp0s(xLeft1, a2.p1.target, a2.b0))
        return (a2New, a1New, TransposeType.LEFT)

    xRight1 = removePrefixMol0(comp0(a2.a0, a2.p1.source), a1.a0)
    xRight2 = removeSuffixMol0(comp0(a1.p1.target, a1.b0), a2.b0)

    if xRight1 is not None and xRight1 == xRight2:
        a2New = fAtom1(a2.a0, a2.p1, comp0s(xRight1, a1.p1.source, a1.b0))
        a1New = fAtom1(comp0s(a2.a0, a2.p1.target, xRight1), a1.p1, a1.b0)
        return (a2New, a1New, TransposeType.RIGHT)

    return (None, None, TransposeType.NONE)

# Returns a pair (newMol1, transposition type), where newMol1 may be None.
def transposeAtom1sAtIndex(mol1, k):
    a1 = mol1.atom1s[k]
    a2 = mol1.atom1s[k+1]
    (newA2, newA1, t) = transposeAtom1s(a1, a2)
    if t == TransposeType.NONE:
        return (None, TransposeType.NONE)
    return (fNonIdMol1(take(mol1.atom1s, k) + (newA2, newA1) + drop(mol1.atom1s, k + 2)), t)

def splitMol1InHalf(m):
    l = len(m.atom1s)
    p = int(l / 2)
    leftAtom1s = m.atom1s[:p]
    rightAtom1s = m.atom1s[p:]
    if not leftAtom1s:
        left = fIdMol1(m.source)
    else:
        left = fNonIdMol1(leftAtom1s)
    if not rightAtom1s:
        right = fIdMol1(m.target)
    else:
        right = fNonIdMol1(rightAtom1s)
    return (left,right)


eqMol1ByMol1 = {}

# Split in half, find all equivalent mol1s on each side, then optionally transpose across the joint.
# Every sequence of transpositions should be found by such, and the split allows us to cache
# sub results.
def stepEquivalentMol1s(mol1):
    if len(mol1) <= 1:
        return set()
    (left, right) = splitMol1InHalf(mol1)

    leftEqs = fEqMol1(left).mol1s
    rightEqs = fEqMol1(right).mol1s
    combined = [comp1(x,y) for x in leftEqs for y in rightEqs]

    stepped = set(combined)
    if len(mol1) >= 2:
        transposeIndex = len(left)
        for c in combined:
            a1 = c.atom1s[transposeIndex - 1]
            a2 = c.atom1s[transposeIndex]
            (a2New, a1New, transposeType) = transposeAtom1s(a1,a2)
            if transposeType != TransposeType.NONE:
                result = fNonIdMol1(take(c.atom1s, transposeIndex - 1) + (a2New, a1New) + drop(c.atom1s, transposeIndex + 1))
                stepped.add(result)
    return stepped


def findEquivalentMol1s(mol1):
    if isinstance(mol1, IdMol1):
        return frozenset([mol1])
    frontier = [mol1]
    equiv = frozenset([mol1])
    while frontier:
        nxt = frontier.pop()
        new = stepEquivalentMol1s(nxt)
        remNew = new.difference(equiv)
        frontier = frontier + list(remNew)
        equiv = equiv.union(remNew)
    return equiv


def fEqMol1(mol1):
    if mol1 in eqMol1ByMol1:
        return eqMol1ByMol1[mol1]
    equiv = findEquivalentMol1s(mol1)
    eqMol1 = EqMol1(equiv)
    for e in equiv:
        eqMol1ByMol1[e] = eqMol1
    return eqMol1


atom2Repo = {}
fAtom2 = buildKeyedBy(atom2Repo, Atom2)
aeIdMol2Repo = {}
fAEIdMol2 = buildKeyedBy(aeIdMol2Repo, AEIdMol2)
aeNonIdMol2Repo = {}
fAENonIdMol2 = buildKeyedBy(aeNonIdMol2Repo, AENonIdMol2)

eqAtom2ByAtom2 = {}

def collapseToAtom2(mol1):
    for ii in range(len(mol1.atom1s)):
        m = mol1.atom1s[ii]
        if isinstance(m.p1, CollapsePrim1):
            return fAtom2(takeMol1(mol1, ii), m.a0, m.p1.prim2, m.b0, dropMol1(mol1, ii+1))


def findEquivalentAtom2s(atom2):
    eqCollapse = fEqMol1(atom2.collapse)
    return frozenset(map(collapseToAtom2, eqCollapse.mol1s))


def ensureMol0(x):
    if isinstance(x, Prim0):
        x = prim0ToMol0(x)
    return x

def ensureEqMol1(x):
    if isinstance(x, Prim1):
        x = prim1ToMol1(x)
    elif isinstance(x, Atom1):
        x = fNonIdMol1((x,))
    if isinstance(x, Mol1):
        x = fEqMol1(x)
    return x

def ensureEqAEMol2(x):
    if isinstance(x, Prim2):
        x = prim2ToAEMol2(x)
    elif isinstance(x, Atom2):
        x = fAENonIdMol2((fEqAtom2(x),))
    elif isinstance(x, EqAtom2):
        x = fAENonIdMol2((x,))
    if isinstance(x, AEMol2):
        x = fEqAEMol2(x)
    return x

def fEqAtom2(atom2):
    if atom2 in eqAtom2ByAtom2:
        return eqAtom2ByAtom2[atom2]
    equiv = findEquivalentAtom2s(atom2)
    eqAtom2 = EqAtom2(equiv)
    for e in equiv:
        eqAtom2ByAtom2[e] = eqAtom2
    return eqAtom2




# If prefix is a prefix of mol1, returns mol1 with the prefix removed.  Otherwise, returns None.
# Does not check equivalent molecules.
def removePrefixMol1(prefix, mol1):
    if isinstance(prefix, IdMol1):
        return mol1

    p = prefix.atom1s
    m = mol1.atom1s
    if m[:len(p)] == p:
        rest = m[len(p):]
        if len(rest) == 0:
            return fIdMol1(mol1.target)
        else:
            return fNonIdMol1(rest)
    else:
        return None

# If suffix is a suffix of mol1, returns mol1 with the suffix removed.  Otherwise, returns None.
# Does not check equivalent molecules.
def removeSuffixMol1(suffix, mol1):
    if isinstance(suffix, IdMol1):
        return mol1

    s = suffix.atom1s
    m = mol1.atom1s
    if m[-len(s):] == s:
        rest = m[:-len(s)]
        if len(rest) == 0:
            return fIdMol1(mol1.source)
        else:
            return fNonIdMol1(rest)
    else:
        return None

# Returns a tuple of (a2',a1') or None
def leftTransposeEqAtom2s(a1,a2):
    a1r = a1.righthand
    a2l = a2.lefthand

    a1rEq = fEqMol1(a1r.r1)
    a2lEq = fEqMol1(a2l.l1)

    # TODO: Move composition outside loops to speed up
    for rightMol in a1rEq.mol1s:
        xLeft1 = removeSuffixMol1(comp1(comp0s(a2l.a0, a2l.p2.source, a2l.b0), a2l.r1), rightMol)
        if xLeft1 is not None:
            for leftMol in a2lEq.mol1s:
                xLeft2 = removePrefixMol1(comp1(a1r.l1, comp0s(a1r.a0, a1r.p2.target, a1r.b0)), leftMol)
                if xLeft1 == xLeft2:
                    a2New = fEqAtom2(fAtom2(comp1s(a1r.l1, comp0s(a1r.a0, a1r.p2.source, a1r.b0), xLeft1), a2l.a0, a2l.p2, a2l.b0, a2l.r1))
                    a1New = fEqAtom2(fAtom2(a1r.l1, a1r.a0, a1r.p2, a1r.b0, comp1s(xLeft1, comp0s(a2l.a0, a2l.p2.target, a2l.b0), a2l.r1)))
                    return (a2New, a1New)
    return (None, None)

# Returns a tuple of (a2',a1') or (None, None)
def rightTransposeEqAtom2s(a1,a2):
    a1l = a1.lefthand
    a2r = a2.righthand

    a1lEq = fEqMol1(a1l.l1)
    a2rEq = fEqMol1(a2r.r1)

    # TODO: Move composition outside loops to speed up
    for leftMol in a1lEq.mol1s:
        xRight1 = removePrefixMol1(comp1(a2r.l1, comp0s(a2r.a0, a2r.p2.source, a2r.b0)), leftMol)
        if xRight1 is not None:
            for rightMol in a2rEq.mol1s:
                xRight2 = removeSuffixMol1(comp1(comp0s(a1l.a0, a1l.p2.target, a1l.b0), a1l.r1), rightMol)
                if xRight1 == xRight2:
                    a2New = fEqAtom2(fAtom2(a2r.l1, a2r.a0, a2r.p2, a2r.b0, comp1s(xRight1, comp0s(a1l.a0, a1l.p2.source, a1l.b0), a1l.r1)))
                    a1New = fEqAtom2(fAtom2(comp1s(a2r.l1, comp0s(a2r.a0, a2r.p2.target, a2r.b0), xRight1), a1l.a0, a1l.p2, a1l.b0, a1l.r1))
                    return (a2New, a1New)
    return (None, None)

def leftTransposeEqAtom2sAtIndex(aeMol2, k):
    a1 = aeMol2.eqAtom2s[k]
    a2 = aeMol2.eqAtom2s[k+1]
    (newA2, newA1) = leftTransposeEqAtom2s(a1, a2)
    if newA2 is None:
        return None
    return fAENonIdMol2(take(aeMol2.eqAtom2s, k) + (newA2, newA1) + drop(aeMol2.eqAtom2s, k + 2))


def splitAEMol2InHalf(m):
    l = len(m.eqAtom2s)
    p = int(l / 2)
    leftAtom2s = m.eqAtom2s[:p]
    rightAtom2s = m.eqAtom2s[p:]
    if not leftAtom2s:
        left = fAEIdMol2(m.source)
    else:
        left = fAENonIdMol2(leftAtom2s)
    if not rightAtom2s:
        right = fAEIdMol2(m.target)
    else:
        right = fAENonIdMol2(rightAtom2s)
    return (left,right)

eqAEMol2ByAEMol2 = {}

def stepEquivalentAEMol2s(aeMol2):
    (left, right) = splitAEMol2InHalf(aeMol2)
    if left == aeMol2 or right == aeMol2:
        return set()

    leftEqs = fEqAEMol2(left).aeMol2s
    rightEqs = fEqAEMol2(right).aeMol2s
    combined = [comp2(x,y) for x in leftEqs for y in rightEqs]

    stepped = set(combined)
    if len(aeMol2) >= 2:
        transposeIndex = len(left)
        for c in combined:
            a1 = c.eqAtom2s[transposeIndex - 1]
            a2 = c.eqAtom2s[transposeIndex]
            (a2NewL, a1NewL) = leftTransposeEqAtom2s(a1,a2)
            if a2NewL:
                result = fAENonIdMol2(take(c.eqAtom2s, transposeIndex - 1) + (a2NewL, a1NewL) + drop(c.eqAtom2s, transposeIndex + 1))
                stepped.add(result)
            (a2NewR, a1NewR) = rightTransposeEqAtom2s(a1,a2)
            if a2NewR:
                result = fAENonIdMol2(take(c.eqAtom2s, transposeIndex - 1) + (a2NewR, a1NewR) + drop(c.eqAtom2s, transposeIndex + 1))
                stepped.add(result)
    return stepped


def findEquivalentAEMol2s(aeMol2):
    if isinstance(aeMol2, AEIdMol2):
        return frozenset([aeMol2])
    frontier = [aeMol2]
    equiv = frozenset([aeMol2])
    while frontier:
        nxt = frontier.pop()
        new = stepEquivalentAEMol2s(nxt)
        remNew = new.difference(equiv)
        frontier = frontier + list(remNew)
        equiv = equiv.union(remNew)
    return equiv


def fEqAEMol2(aeMol2):
    if aeMol2 in eqAEMol2ByAEMol2:
        return eqAEMol2ByAEMol2[aeMol2]
    equiv = findEquivalentAEMol2s(aeMol2)
    eqAEMol2 = EqAEMol2(equiv)
    for e in equiv:
        eqAEMol2ByAEMol2[e] = eqAEMol2
    return eqAEMol2

def primSetEqMol1(eqMol1):
    inst = next(iter(eqMol1.mol1s))
    if isinstance(inst, IdMol1):
        return set()
    else:
        pset = set()
        for atom in inst.atom1s:
            pset.add(atom.p1)
        return pset
