from ontology import *

# Gives a set of non-trivial decompositions of mol0 via tensor, as pairs (m1, m2)
@memoize1
def tensorDecompMol0(mol0):
    decomps = set()
    # If len(mol0) <= 1, this does nothing.
    for ii in range(1, len(mol0)):
        left = fMol0(mol0.prim0s[:ii])
        right = fMol0(mol0.prim0s[ii:])
        decomps.add((left, right))
    return frozenset(decomps)

###########################################################
# Horizontal decomposition, dimension 1

# Returns a list of all non-trivial horizontal decompositions of the eqMol1, as a set of tuples.
@memoize1
def horizontalDecompEqMol1(eqMol1):
    decomps = set()
    for mol1 in eqMol1.mol1s:
        for ii in range(1, len(mol1)):
            decomps.add((fEqMol1(takeMol1(mol1, ii)), fEqMol1(dropMol1(mol1, ii))))
    return frozenset(decomps)


###########################################################
# Tensor decomposition, dimension 1

# Returns a list of all non-trivial prefixes of a tuple (including the whole tuple).
def nonTrivPrefixes(t):
    prefixes = []
    for ii in range(1,len(t) + 1):
        prefixes.append(take(t, ii))
    return prefixes

def nonTrivSuffixes(t):
    suffixes = []
    for ii in range(0,len(t)):
        suffixes.append(drop(t, ii))
    return suffixes

# Gives a set of non-trivial decompositions of eqMol1 via tensor, as pairs (m1, m2).
@memoize1
def tensorDecompEqMol1(eqMol1):
    decomps = set()
    if len(eqMol1) == 0:
        # Identity decomps
        mol0Decomps = tensorDecompMol0(eqMol1.source)
        for (left, right) in mol0Decomps:
            decomps.add((fEqMol1(fIdMol1(left)), fEqMol1(fIdMol1(right))))
    else:
        # Prefix decomps
        # Shared prefixes are independent of transpositions, so we only need to check one instance.
        mol1Inst = next(iter(eqMol1.mol1s))
        decomps = decomps.union(map(lambda x: (fEqMol1(x[0]), fEqMol1(x[1])), idPrefixDecompMol1(mol1Inst)))

        # Suffix decomps
        decomps = decomps.union(map(lambda x: (fEqMol1(x[0]), fEqMol1(x[1])), idSuffixDecompMol1(mol1Inst)))

        # Block decomps
        decomps = decomps.union(map(lambda x: (fEqMol1(x[0]), fEqMol1(x[1])),blockTensorDecompEqMol1(eqMol1)))
    return frozenset(decomps)

# Gives a set of non-trivial decompositions of a NonIdMol1 via tensor, where the first of the pair is an identity 1-molecule
def idPrefixDecompMol1(mol1):
    decomps = set()
    # Potential prefixes are the prefixes of the first atom.
    firstAtom = mol1.atom1s[0]
    prefixMol0s = map(lambda x: fMol0(x), nonTrivPrefixes(firstAtom.a0.prim0s))
    for pre in prefixMol0s:
        # Simultaneously check if it's a prefix and build version of the atoms without it.
        newAtoms = []
        broken = False
        for atom in mol1.atom1s:
            removed = removePrefixMol0(pre, atom.a0)
            if removed is not None:
                newAtoms.append(fAtom1(removed, atom.p1, atom.b0))
            else:
                broken = True
                break
        if not(broken):
            decomps.add((fIdMol1(pre), fNonIdMol1(tuple(newAtoms))))
    return decomps

# Gives a set of non-trivial decompositions of a NonIdMol1 via tensor, where the second of the pair is an identity 1-molecule
def idSuffixDecompMol1(mol1):
    decomps = set()
    firstAtom = mol1.atom1s[0]
    suffixMol0s = map(lambda x: fMol0(x), nonTrivSuffixes(firstAtom.b0.prim0s))
    for suf in suffixMol0s:
        newAtoms = []
        broken = False
        for atom in mol1.atom1s:
            removed = removeSuffixMol0(suf, atom.b0)
            if removed is not None:
                newAtoms.append(fAtom1(atom.a0, atom.p1, removed))
            else:
                broken = True
                break
        if not(broken):
            decomps.add((fNonIdMol1(tuple(newAtoms)), fIdMol1(suf)))
    return decomps

def blockTensorDecompEqMol1(eqMol1):
    decomps = set()
    for mol1 in eqMol1.mol1s:
        for ii in range(1,len(mol1)):
            decomps = decomps.union(blockTensorDecompMol1(mol1, ii))
    return decomps

# Returns a pair (m1, m2) with m2 a maximal identity molecule (perhaps 1_\emptymol)
def maxSuffixDecompMol1(mol1):
    suffixDecomps = idSuffixDecompMol1(mol1)
    if len(suffixDecomps) > 0:
        maximal = max(suffixDecomps, key=lambda x: len(x[1].mol0))
        m1 = maximal[0]
        m2 = maximal[1]
    else:
        m1 = mol1
        m2 = fIdMol1(fMol0(()))
    return (m1, m2)

# Similar to above
def maxPrefixDecompMol1(mol1):
    prefixDecomps = idPrefixDecompMol1(mol1)
    if len(prefixDecomps) > 0:
        maximal = max(prefixDecomps, key=lambda x: len(x[0].mol0))
        m1 = maximal[0]
        m2 = maximal[1]
    else:
        m1 = fIdMol1(fMol0(()))
        m2 = mol1
    return (m1, m2)

# Attempts to decompose the given molecule1 as a block tensor at index k.  Returns None if fails.  If succeeds, returns a pair (m1, m2).
def blockTensorDecompMol1(mol1, k):
    m = mol1
    for jj in range(0, len(mol1) - k):
        for ii in range(1,k+1):
            (m, t) = transposeAtom1sAtIndex(m, k + jj - ii)
            if t == TransposeType.RIGHT or t == TransposeType.NONE:
                return set()

    # If that worked, we're of the right tensor form.
    # We can read the tensor form off by looking at the maximal unchanged blocks between the interchanges.
    block1Pre = fNonIdMol1(take(mol1.atom1s, k))
    block2Pre = fNonIdMol1(drop(mol1.atom1s, k))

    (g, maxSuffix) = maxSuffixDecompMol1(block1Pre)
    (discard, h) = maxPrefixDecompMol1(block2Pre)
    w = removeSuffixMol0(h.source, maxSuffix.mol0)

    decomps = set()
    for ii in range(0,len(w) + 1):
        left = comp0(g, takeMol0(w, ii))
        right = comp0(dropMol0(w, ii), h)
        decomps.add((left, right))
    return decomps


###########################################################
# Vertical decomposition, dimension 2

# Returns a list of all non-trivial horizontal decompositions of the eqMol1, as a set of tuples.
@memoize1
def verticalDecompEqAEMol2(eqAEMol2):
    decomps = set()
    for aeMol2 in eqAEMol2.aeMol2s:
        for ii in range(1, len(aeMol2)):
            decomps.add((fEqAEMol2(takeAEMol2(aeMol2, ii)), fEqAEMol2(dropAEMol2(aeMol2, ii))))
    return frozenset(decomps)


##################################################
# Horizontal decomposition, dimension 2

# Returns a pair (m1, m2) with m2 a maximal identity molecule (perhaps 1_{1_cod})
def maxSuffixDecompAEMol2(aeMol2):
    suffixDecomps = idSuffixDecompAEMol2(aeMol2)
    if len(suffixDecomps) > 0:
        maximal = max(suffixDecomps, key=lambda x: len(x[1].eqMol1))
        m1 = maximal[0]
        m2 = maximal[1]
    else:
        m1 = aeMol2
        m2 = fAEIdMol2(fIdMol1(aeMol2.target.target))
    return (m1, m2)

# Similar to above
def maxPrefixDecompAEMol2(aeMol2):
    prefixDecomps = idPrefixDecompAEMol2(aeMol2)
    if len(prefixDecomps) > 0:
        maximal = max(prefixDecomps, key=lambda x: len(x[0].eqMol1))
        m1 = maximal[0]
        m2 = maximal[1]
    else:
        m1 = fAEIdMol2(fIdMol1(aeMol2.source.source))
        m2 = aeMol2
    return (m1, m2)


def idPrefixDecompAEMol2(aeMol2):
    decomps = set()

    prefixes = []
    for atom in aeMol2.eqAtom2s[0].atom2s:
        prefixes += map(fNonIdMol1, nonTrivPrefixes(atom.l1.atom1s))

    # Enough to cancel the prefix in one equivalent atom, by cancellation and collapse definition of equality.
    for pre in prefixes:
        newEqAtoms = []
        broke = False
        for ii in range(len(aeMol2)):
            found = False
            for atom in aeMol2.eqAtom2s[ii].atom2s:
                newL1 = removePrefixMol1(pre, atom.l1)
                if newL1 is not None:
                    found = True
                    newEqAtoms.append(fEqAtom2(fAtom2(newL1, atom.a0, atom.p2, atom.b0, atom.r1)))
                    break
            if not(found):
                broke = True
                break
        if broke:
            continue
        decomps.add((fAEIdMol2(fEqMol1(pre)), fAENonIdMol2(tuple(newEqAtoms))))
    return decomps

def idSuffixDecompAEMol2(aeMol2):
    decomps = set()

    suffixes = []
    for atom in aeMol2.eqAtom2s[0].atom2s:
        suffixes += map(fNonIdMol1, nonTrivSuffixes(atom.r1.atom1s))

    # Enough to cancel the suffix in one equivalent atom, by cancellation and collapse definition of equality.
    for suf in suffixes:
        newEqAtoms = []
        broke = False
        for ii in range(len(aeMol2)):
            found = False
            for atom in aeMol2.eqAtom2s[ii].atom2s:
                newR1 = removeSuffixMol1(suf, atom.r1)
                if newR1 is not None:
                    found = True
                    newEqAtoms.append(fEqAtom2(fAtom2(atom.l1, atom.a0, atom.p2, atom.b0, newR1)))
                    break
            if not(found):
                broke = True
                break
        if broke:
            continue
        decomps.add((fAENonIdMol2(tuple(newEqAtoms)), fAEIdMol2(fEqMol1(suf))))
    return decomps

def blockHorizontalDecompEqAEMol2(m):
    decomps = set()
    for aeMol2 in m.aeMol2s:
        for ii in range(1,len(m)):
            decomps = decomps.union(blockHorizontalDecompAEMol2(aeMol2, ii))
    return decomps

def blockHorizontalDecompAEMol2(aeMol2, k):
    m = aeMol2
    for jj in range(0, len(aeMol2) - k):
        for ii in range(1,k+1):
            m = leftTransposeEqAtom2sAtIndex(m, k + jj - ii)
            if m is None:
                return set()

    block1Pre = fAENonIdMol2(take(aeMol2.eqAtom2s, k))
    block2Pre = fAENonIdMol2(drop(aeMol2.eqAtom2s, k))

    (g, maxSuffix) = maxSuffixDecompAEMol2(block1Pre)
    (discard, h) = maxPrefixDecompAEMol2(block2Pre)

    # Need to remove the source of h as a suffix on maxSuffix, where we may need to change up to equivalence.
    hSourceEqMol1 = h.source
    suffixEqMol1 = maxSuffix.eqMol1
    w = None
    for x in suffixEqMol1.mol1s:
        found = False
        for y in hSourceEqMol1.mol1s:
            w = removeSuffixMol1(y, x)
            if w is not None:
                found = True
                break
        if found:
            break

    wEqMol1 = fEqMol1(w)
    wPartitions = set()
    for wMol1 in wEqMol1.mol1s:
        for ii in range(0,len(wMol1) + 1):
            wPartitions.add((takeMol1(wMol1, ii), dropMol1(wMol1, ii)))

    decomps = set()
    for (wLeft, wRight) in wPartitions:
        left = comp1(g, wLeft)
        right = comp1(wRight, h)
        decomps.add((left, right))
    return decomps

# Gives a set of non-trivial decompositions of eqAEMol2 via horizontal composition, as pairs (m1, m2).
@memoize1
def horizontalDecompEqAEMol2(eqAEMol2):
    decomps = set()
    if len(eqAEMol2) == 0:
        # Identity decomps
        mol1Decomps = horizontalDecompEqMol1(eqAEMol2.source)
        for (left, right) in mol1Decomps:
            decomps.add((fEqAEMol2(fAEIdMol2(left)), fEqAEMol2(fAEIdMol2(right))))
    else:
        # TODO: Haven't verified that molecule2 1-prefixes are independent of transpositions.
        # For now, just look over all instances.
        for aeMol2Inst in eqAEMol2.aeMol2s:
            decomps = decomps.union(map(lambda x: (fEqAEMol2(x[0]), fEqAEMol2(x[1])), idPrefixDecompAEMol2(aeMol2Inst)))
            decomps = decomps.union(map(lambda x: (fEqAEMol2(x[0]), fEqAEMol2(x[1])), idSuffixDecompAEMol2(aeMol2Inst)))

        decomps = decomps.union(map(lambda x: (fEqAEMol2(x[0]), fEqAEMol2(x[1])),blockHorizontalDecompEqAEMol2(eqAEMol2)))
    return frozenset(decomps)

##################################################
# Tensor decomposition, dimension 2

def isLeftTensorTransposableAtIndex(aeMol2, k):
    return areLeftTensorTransposable(aeMol2.eqAtom2s[k], aeMol2.eqAtom2s[k+1])

# Returns true if the two are left tensor transposable.
def areLeftTensorTransposable(eqA1, eqA2):
    a1rBase = eqA1.righthand
    a2lBase = eqA2.lefthand

    # Temporarily replace with with distinct primitives; otherwise this fails when the primitive of both atoms is the same,
    # since they can't be distinguished in the shared collapse.
    p1 = ConstPrim2("", a1rBase.p2.source, a1rBase.p2.target)
    p2 = ConstPrim2("", a2lBase.p2.source, a2lBase.p2.target)
    a1r = fAtom2(a1rBase.l1, a1rBase.a0, p1, a1rBase.b0, a1rBase.r1)
    a2l = fAtom2(a2lBase.l1, a2lBase.a0, p2, a2lBase.b0, a2lBase.r1)

    a1rEq = fEqMol1(a1r.r1)
    a2lEq = fEqMol1(a2l.l1)

    # TODO: Copied from left tensor transpose, factor out of both.
    for rightMol in a1rEq.mol1s:
        xLeft1 = removeSuffixMol1(comp1(comp0s(a2l.a0, a2l.p2.source, a2l.b0), a2l.r1), rightMol)
        if xLeft1 is not None:
            for leftMol in a2lEq.mol1s:
                xLeft2 = removePrefixMol1(comp1(a1r.l1, comp0s(a1r.a0, a1r.p2.target, a1r.b0)), leftMol)
                if xLeft1 == xLeft2:
                    sharedCollapse = comp1s(a1r.l1, comp0s(a1r.a0, a1r.p2.collapse, a1r.b0), xLeft1, comp0s(a2l.a0, a2l.p2.collapse, a2l.b0), a2l.r1)
                    tensorDecomps = tensorDecompEqMol1(fEqMol1(sharedCollapse))
                    validDecomps = []
                    for (left, right) in tensorDecomps:
                        if a1r.p2.collapse in primSetEqMol1(left) and a2l.p2.collapse in primSetEqMol1(right):
                            #validDecomps.append((left, right))
                            return True
                    #(bestLeft, bestRight) = max(validDecomps, key=lambda x: len(source0(x[0])))
                    #return (fEqAtom2(collapseToAtom2(bestLeft)), fEqAtom2(collapseToAtom2(bestRight)))
    return False

def nonTrivTensorPrefixes(eqAtom2):
    atom2 = next(iter(eqAtom2.atom2s))
    decomps = tensorDecompEqMol1(fEqMol1(atom2.collapse))
    prefixes = []
    for (left, right) in decomps:
        if atom2.p2.collapse in primSetEqMol1(right):
            prefixes.append(left)
    return prefixes

def nonTrivTensorSuffixes(eqAtom2):
    atom2 = next(iter(eqAtom2.atom2s))
    decomps = tensorDecompEqMol1(fEqMol1(atom2.collapse))
    suffixes = []
    for (left, right) in decomps:
        if atom2.p2.collapse in primSetEqMol1(left):
            suffixes.append(right)
    return suffixes

def removeTensorPrefixEqMol1(preEqMol1, eqMol1):
    if len(preEqMol1) == 0 and preEqMol1.source == fMol0(()):
        return eqMol1
    if preEqMol1 == eqMol1:
        return fEqMol1(fIdMol1(fMol0(())))

    decomps = tensorDecompEqMol1(eqMol1)
    for (left, right) in decomps:
        if left == preEqMol1:
            return right
    return None

def removeTensorSuffixEqMol1(sufEqMol1, eqMol1):
    if len(sufEqMol1) == 0 and sufEqMol1.source == fMol0(()):
        return eqMol1
    if sufEqMol1 == eqMol1:
        return fEqMol1(fIdMol1(fMol0(())))

    decomps = tensorDecompEqMol1(eqMol1)
    for (left, right) in decomps:
        if right == sufEqMol1:
            return left
    return None

def removeTensorPrefixEqAtom2(preEqMol1, eqAtom2):
    if len(preEqMol1) == 0 and preEqMol1.source == fMol0(()):
        return eqAtom2

    decomps = tensorDecompEqMol1(fEqMol1(next(iter(eqAtom2.atom2s)).collapse))
    for (left, right) in decomps:
        if left == preEqMol1:
            return fEqAtom2(collapseToAtom2(next(iter(right.mol1s))))
    return None

def removeTensorSuffixEqAtom2(sufEqMol1, eqAtom2):
    if len(sufEqMol1) == 0 and sufEqMol1.source == fMol0(()):
        return eqAtom2

    decomps = tensorDecompEqMol1(fEqMol1(next(iter(eqAtom2.atom2s)).collapse))
    for (left, right) in decomps:
        if right == sufEqMol1:
            return fEqAtom2(collapseToAtom2(next(iter(left.mol1s))))
    return None

def idTensorPrefixDecompAEMol2(aeMol2):
    decomps = set()

    # Tensor prefixes are independent of choices of atom2 instance.
    prefixes = nonTrivTensorPrefixes(aeMol2.eqAtom2s[0])
    for pre in prefixes:
        newEqAtoms = []
        broke = False
        for eqAtom in aeMol2.eqAtom2s:
            result = removeTensorPrefixEqAtom2(pre, eqAtom)
            if result is None:
                broke = True
                break
            newEqAtoms.append(result)
        if not(broke):
            decomps.add((fAEIdMol2(pre), fAENonIdMol2(tuple(newEqAtoms))))
    return decomps

def idTensorSuffixDecompAEMol2(aeMol2):
    decomps = set()

    # Tensor suffixes are independent of choices of atom2 instance.
    suffixes = nonTrivTensorSuffixes(aeMol2.eqAtom2s[0])
    for suf in suffixes:
        newEqAtoms = []
        broke = False
        for eqAtom in aeMol2.eqAtom2s:
            result = removeTensorSuffixEqAtom2(suf, eqAtom)
            if result is None:
                broke = True
                break
            newEqAtoms.append(result)
        if not(broke):
            decomps.add((fAENonIdMol2(tuple(newEqAtoms)), fAEIdMol2(suf)))
    return decomps

# Returns a pair (m1, m2) with m1 a maximal identity molecule (perhaps 1_{1_{1_emptymol}})
def maxTensorPrefixDecompAEMol2(aeMol2):
    prefixDecomps = idTensorPrefixDecompAEMol2(aeMol2)
    if len(prefixDecomps) > 0:
        maximal = max(prefixDecomps, key=lambda x: len(source0(x[0])))
        m1 = maximal[0]
        m2 = maximal[1]
    else:
        m1 = fAEIdMol2(ensureEqMol1(fMol0(())))
        m2 = aeMol2
    return (m1, m2)

# Returns a pair (m1, m2) with m2 a maximal identity molecule (perhaps 1_{1_{1_emptymol}})
def maxTensorSuffixDecompAEMol2(aeMol2):
    suffixDecomps = idTensorSuffixDecompAEMol2(aeMol2)
    if len(suffixDecomps) > 0:
        maximal = max(suffixDecomps, key=lambda x: len(source0(x[1])))
        m1 = maximal[0]
        m2 = maximal[1]
    else:
        m1 = aeMol2
        m2 = fAEIdMol2(ensureEqMol1(fMol0(())))
    return (m1, m2)

def blockTensorDecompEqAEMol2(m):
    decomps = set()
    for aeMol2 in m.aeMol2s:
        for ii in range(1,len(m)):
            decomps = decomps.union(blockTensorDecompAEMol2(aeMol2, ii))
    return decomps

def blockTensorDecompAEMol2(aeMol2, k):
    m = aeMol2
    for jj in range(0, len(aeMol2) - k):
        for ii in range(1,k+1):
            if isLeftTensorTransposableAtIndex(m, k + jj - ii):
                m = leftTransposeEqAtom2sAtIndex(m, k + jj - ii)
            else:
                return set()

    block1Pre = fAENonIdMol2(take(aeMol2.eqAtom2s, k))
    block2Pre = fAENonIdMol2(drop(aeMol2.eqAtom2s, k))

    (g, maxSuffix) = maxTensorSuffixDecompAEMol2(block1Pre)
    (discard, h) = maxTensorPrefixDecompAEMol2(block2Pre)

    w = removeTensorSuffixEqMol1(h.source, maxSuffix.eqMol1)
    wDecomps = set(tensorDecompEqMol1(w))
    wDecomps.add((w, fEqMol1(fIdMol1(fMol0(())))))
    wDecomps.add((fEqMol1(fIdMol1(fMol0(()))), w))

    decomps = set()
    for (wLeft, wRight) in wDecomps:
        left = comp0(g, wLeft)
        right = comp0(wRight, h)
        decomps.add((left, right))
    return decomps

# Gives a set of non-trivial decompositions of eqAEMol2 via tensor, as pairs (m1, m2).
@memoize1
def tensorDecompEqAEMol2(eqAEMol2):
    decomps = set()
    if len(eqAEMol2) == 0:
        # Identity decomps
        mol1Decomps = tensorDecompEqMol1(eqAEMol2.source)
        for (left, right) in mol1Decomps:
            decomps.add((fEqAEMol2(fAEIdMol2(left)), fEqAEMol2(fAEIdMol2(right))))
    else:
        # TODO: Haven't verified that molecule2 1-prefixes are independent of transpositions.
        # For now, just look over all instances.
        for aeMol2Inst in eqAEMol2.aeMol2s:
            decomps = decomps.union(map(lambda x: (fEqAEMol2(x[0]), fEqAEMol2(x[1])), idTensorPrefixDecompAEMol2(aeMol2Inst)))
            decomps = decomps.union(map(lambda x: (fEqAEMol2(x[0]), fEqAEMol2(x[1])), idTensorSuffixDecompAEMol2(aeMol2Inst)))

        decomps = decomps.union(map(lambda x: (fEqAEMol2(x[0]), fEqAEMol2(x[1])),blockTensorDecompEqAEMol2(eqAEMol2)))
    return frozenset(decomps)
