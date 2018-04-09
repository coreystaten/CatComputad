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

    def __repr__(self):
        return "Mol0(" + ", ".join(map(repr, self.prim0s)) + ")"

    def __len__(self):
        return len(self.prim0s)


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
        self.name = "collapse(" + prim2.name + ")"
        self.source = prim2.source.source
        self.target = prim2.target.target

    def __repr__(self):
        return self.name

class Atom1(IdHashed):
    def __init__(self,a0,p1,b0):
        self.a0 = a0
        self.b0 = b0
        self.p1 = p1
        self.source = comp0s([self.a0, self.p1.source, self.b0])
        self.target = comp0s([self.a0, self.p1.target, self.b0])

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

    def __repr__(self):
        return "IdMol1(" + repr(self.mol0) + ")"

    def __len__(self):
        return 0


class NonIdMol1(Mol1):
    def __init__(self, atom1s):
        self.atom1s = atom1s
        self.source = self.atom1s[0].source
        self.target = self.atom1s[-1].target

    def __repr__(self):
        return "Mol1(" + ", ".join(map(repr, self.atom1s)) + ")"

    def __len__(self):
        return len(self.atom1s)

class EqMol1(IdHashed):
    def __init__(self, mol1s):
        self.mol1s = mol1s
        self.length = len(list(self.mol1s)[0])

    def __repr__(self):
        return "EqMol1(" + repr(list(self.mol1s)) + ")"

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

    def __repr__(self):
        return self.name

class Atom2(IdHashed):
    def __init__(self,l1,a0,p2,b0,r1):
        self.l1 = l1
        self.a0 = a0
        self.p2 = p2
        self.b0 = b0
        self.r1 = r1
        self.source = comp1s([self.l1, comp0s([self.a0, self.p2.source, self.b0]), self.r1])
        self.target = comp1s([self.l1, comp0s([self.a0, self.p2.target, self.b0]), self.r1])
        self.collapse = comp1s([self.l1, comp0s([self.a0, self.p2.collapse, self.b0]), self.r1])

    def __repr__(self):
        return "Atom2(" + repr(self.l1) + ", " + repr(self.a0) + ", " + repr(self.p1) + ", " + repr(self.b0) + ", " + repr(self.r1) + ")"

class EqAtom2(IdHashed):
    def __init__(self,atom2s):
        self.atom2s = atom2s
        self.righthand = min(self.atom2s, key=lambda x: len(x.l1.atom1s))
        self.lefthand = min(self.atom2s, key=lambda x: len(x.r1.atom1s))
        self.source = self.atom2s[0].source
        self.target = self.atom2s[0].target

    def __repr__(self):
        return "EqAtom2(" + repr(self.atom2s) + ")"

class AEMol2(IdHashed):
    pass

class AEIdMol2(AEMol2):
    def __init__(self, eqMol1):
        self.eqMol1 = eqMol1
        self.atom2s = ()
        self.source = self.eqMol1
        self.target = self.eqMol1

    def __repr__(self):
        return "IdMol2(" + repr(self.eqMol1) + ")"

    def __len__(self):
        return 0

class AENonIdMol2(AEMol2):
    def __init__(self, eqAtom2s):
        self.eqAtom2s = eqAtom2s
        self.source = self.eqAtom2s[0].source
        self.target = self.eqAtom2s[-1].target

    def __repr__(self):
        return "Mol2(" + ", ".join(map(repr, self.eqAtom2s)) + ")"

    def __len__(self):
        return len(self.eqAtom2s)

class EqAEMol2(IdHashed):
    def __init__(self, mol2s):
        self.mol2s = mol2s
        self.length = len(list(self.mol2s)[0])

    def __repr__(self):
        return "EqMol2(" + repr(list(self.mol2s)) + ")"

    def __len__(self):
        return self.length

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

def takeMol1(mol1, n):
    if n <= 0:
        return fIdMol1(mol1.source)
    return fNonIdMol1(take(mol1.atom1s, n))

def dropMol1(mol1, n):
    if n >= len(mol1):
        return fIdMol1(mol1.target)
    return fNonIdMol1(drop(mol1.atom1s, n))

def takeMol2(mol2, n):
    if n <= 0:
        return fIdMol2(mol2.source)
    return fNonIdMol2(take(mol2.atom2s, n))

def dropMol2(mol2, n):
    if n >= len(mol2):
        return fIdMol2(mol2.target)
    return fNonIdMol2(drop(mol2.atom2s, n))

def prim0ToMol0(p):
    return fMol0((p,))

def prim1ToAtom1(p):
    return fAtom1(fMol0(()), p, fMol0(()))

def prim1ToMol1(p):
    return fNonIdMol1((prim1ToAtom1(p),))

def atom1ToMol1(a):
    return fNonIdMol1((a,))


def comp0(x,y):
    if isinstance(x, Prim0):
        x = prim0ToMol0(x)
    if isinstance(y, Prim0):
        y = prim0ToMol0(y)

    if isinstance(x, Mol0) and isinstance(y, Mol0):
        return fMol0(x.prim0s + y.prim0s)

    if isinstance(x, Mol0) and isinstance(y, Atom1):
        return fAtom1(comp0(x, y.a0), y.p1, y.b0)
    if isinstance(x, Atom1) and isinstance(y, Mol0):
        return fAtom1(x.a0, x.p1, comp0(x.b0, y))

    if isinstance(x, Mol0) and isinstance(y, Mol1):
        return fMol1(tuple(map(lambda z: comp0(x, z), y.atom1s)))
    if isinstance(x, Mol1) and isinstance(y, Mol0):
        return fMol1(tuple(map(lambda z: comp0(z, y), x.atom1s)))

    if isinstance(x, Prim1):
        x = prim1ToMol1(x)
    elif isinstance(x, Atom1):
        x = atom1ToMol1(x)
    if isinstance(y, Prim1):
        y = prim1ToMol1(y)
    elif isinstance(y, Atom1):
        y = atom1ToMol1(y)

    if isinstance(x, IdMol1) and isinstance(y, IdMol1):
        return fIdMol1(comp0(x.mol0, y.mol0))

    if isinstance(x, Mol1) and isinstance(y, Mol1):
        firstAtoms = tuple(map(lambda z: comp0(z, y.source), x.atom1s))
        secondAtoms = tuple(map(lambda z: comp0(x.target, z), y.atom1s))
        return fNonIdMol1(firstAtoms + secondAtoms)

    raise Exception("Invalid 0-composition between %s and %s" % (str(type(x)), str(type(y))))

def comp0s(xs):
    result = xs[0]
    for x in xs[1:]:
        result = comp0(result,x)
    return result

def comp1(x,y):
    if isinstance(x, Prim1):
        x = prim1ToMol1(x)
    if isinstance(y, Prim1):
        y = prim1ToMol1(y)
    if isinstance(x, Atom1):
        x = atom1ToMol1(x)
    if isinstance(y, Atom1):
        y = atom1ToMol1(y)

    if isinstance(x, Mol1) and isinstance(y, Mol1):
        if x.target != y.source:
            raise Exception("Source target mismatch in 1-composition: %s and %s" % (repr(x), repr(y)))
        return fNonIdMol1(x.atom1s + y.atom1s)

    raise Exception("Invalid 1-composition between %s and %s" % (str(type(x)), str(type(y))))


def comp1s(xs):
    result = xs[0]
    for x in xs[1:]:
        result = comp1(res,x)
    return result


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
def transposeAtom1s(a1,a2):
    xLeft1 = removeSuffixMol0(comp0(a2.p1.source, a2.b0), a1.b0)
    xLeft2 = removePrefixMol0(comp0(a1.a0, a1.p1.target), a2.a0)

    if xLeft1 is not None and xLeft1 == xLeft2:
        a2New = fAtom1(comp0s([a1.a0, a1.p1.source, xLeft1]), a2.p1, a2.b0)
        a1New = fAtom1(a1.a0, a1.p1, comp0s([xLeft1, a2.p1.target, a2.b0]))
        return (a2New, a1New, TransposeType.LEFT)

    xRight1 = removePrefixMol0(comp0(a2.a0, a2.p1.source), a1.a0)
    xRight2 = removeSuffixMol0(comp0(a1.p1.target, a1.b0), a2.b0)

    if xRight1 is not None and xRight1 == xRight2:
        a2New = fAtom1(a2.a0, a2.p1, comp0s([xRight1, a1.p1.source, a1.b0]))
        a1New = fAtom1(comp0s([a2.a0, a2.p1.target, xRight1]), a1.p1, a1.b0)
        return (a2New, a1New, TransposeType.RIGHT)

    return (None, None, TransposeType.NONE)

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

def stepEquivalentMol1s(mol1):
    (left, right) = splitMol1InHalf(mol1)
    if left == mol1 or right == mol1:
        return set()

    leftEqs = fEqMol1(left).mol1s
    rightEqs = fEqMol1(right).mol1s
    combined = [comp1(x,y) for x in leftEqs for y in rightEqs]

    stepped = set(combined)
    if len(mol1.atom1s) >= 2:
        transposeIndex = len(left.atom1s)
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
            return fAtom2(takeMol1(mol1, ii), comp0s([m.a0, m.p1.prim2, m.b0]), dropMol1(mol1, ii+1))


def findEquivalentAtom2s(atom2):
    eqCollapse = fEqMol1(atom2.collapse)
    return frozenset(map(collapseToAtom2, eqCollapse.mol1s))


def fEqAtom2(atom2):
    if atom2 in eqAtom2ByAtom21:
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
        xLeft1 = removeSuffixMol1(comp1(comp0s([a2l.a0, a2l.p2.source, a2l.b0]), a2l.r1), rightMol)
        if xLeft1:
            for leftMol in a2lEq.mol1s:
                xLeft2 = removePrefixMol1(comp1(a1r.l1, comp0s([a1r.a0, a1r.p2.target, a1r.b0])), leftMol)
                if xLeft1 == xLeft2:
                    a2New = fEqAtom2(fAtom2(comp1s([a1r.l1, comp0s([a1r.a0, a1r.p2.source, a1r.b0]), xLeft1]), a2l.a0, a2l.p2, a2l.b0, a2l.r1))
                    a1New = fEqAtom2(fAtom2(a1r.l1, a1r.a0, a1r.p2, a1r.b0, comp1s([xLeft1, comp0s([a2l.a0, a2l.p2.target, a2l.b0]), a2l.r1))
                    return (a2New, a1New)
    return None

# Returns a tuple of (a2',a1') or None
def rightTransposeEqAtom2s(a1,a2):
    a1l = a1.lefthand
    a2r = a2.righthand

    a1lEq = fEqMol1(a1l.l1)
    a2rEq = fEqMol(a2r.r1)

    # TODO: Move composition outside loops to speed up
    for leftMol in a1lEq.mol1s:
        xRight1 = removePrefixMol1(comp1(a2r.l1, comp0s([a2r.a0, a2r.p2.source, a2r.b0])), leftMol)
        if xRight1:
            for rightMol in a2rEq.mol1s:
                xRight2 = removeSuffixMol1(comp1(comp0s([a1l.a0, a1l.p2.target, a1l.b0]), a1l.r1), rightMol)
                if xRight1 == xRight2:
                    a2New = fEqAtom2(fAtom2(a2r.l1, a2r.a0, a2r.p2, a2r.b0, comp1s([xRight1, comp0s([a1l.a0, a1l.p2.source, a1l.b0]), a1l.r1])))
                    a1New = fEqAtom2(fAtom2(comp1s([a2r.l1, comp0s([a2r.a0, a2r.p2.target, a2r.b0]), xRight1]), a1l.a0, a1l.p2, a1l.b0, a1l.r1))
                    return (a2New, a1New)
    return None

a = ConstPrim0("a")
b = ConstPrim0("b")
c = ConstPrim0("c")
d = ConstPrim0("d")

f = ConstPrim1("f",a,b)
g = ConstPrim1("g",c,d)

l = comp0(f,g)
r = fEqMol1(l)
