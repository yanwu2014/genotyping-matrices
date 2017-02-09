from collections import defaultdict
import re

_DEBUG = 0


# "Fast Text Searching with Errors" by Sun Wu and Udi Manber
# TR 91-11, Dept of Computer Science, University of Arizona, June 1991.
# http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.20.8854

def approxHammingSearch(pattern, text):
    """Generate (Hamming_dist, start_offset)
    for matches with distance 0 or 1"""
    m = len(pattern)
    S_table = defaultdict(int)
    for i, c in enumerate(pattern):
        S_table[c] |= 1 << i
    R0 = 0
    R1 = 0
    mask = 1 << (m - 1)
    for j, c in enumerate(text):
        S = S_table[c]
        shR0 = (R0 << 1) | 1
        R0 = shR0 & S
        R1 = ((R1 << 1) | 1) & S | shR0
        if _DEBUG:
            print "j= %2d msk=%s S=%s R0=%s R1=%s" \
                % tuple([j] + map(bitstr, [mask, S, R0, R1]))
        if R0 & mask: # exact match
            yield 0, j - m + 1
        elif R1 & mask: # match with one substitution
            yield 1, j - m + 1

if _DEBUG:

    def bitstr(num, mlen=8):
       wstr = ""
       for i in xrange(mlen):
          if num & 1:
             wstr = "1" + wstr
          else:
             wstr = "0" + wstr
          num >>= 1
       return wstr

def Ham_dist(s1, s2):
    """Calculate Hamming distance between 2 sequences."""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def long_check(pattern, text):
    """Naively and understandably generate (Hamming_dist, start_offset)
    for matches with distance 0 or 1"""
    m = len(pattern)
    for i in xrange(len(text) - m + 1):
        d = Ham_dist(pattern, text[i:i+m])
        if d < 2:
            yield d, i

def Paul_McGuire_regex(pattern, text):
    searchSeqREStr = (
        '('
        + pattern
        + ')|('
        + ')|('.join(
            pattern[:i]
            + "[ACTGN]".replace(c,'')
            + pattern[i+1:]
            for i,c in enumerate(pattern)
            )
        + ')'
        )
    searchSeqRE = re.compile(searchSeqREStr)
    for match in searchSeqRE.finditer(text):
        locn = match.start()
        dist = int(bool(match.lastindex - 1))
        yield dist, locn


if __name__ == "__main__":

    genome1 = "TTTACGTAAACTAAACTGTAA"
    #         01234567890123456789012345
    #                   1         2

    tests = [
        (genome1, "ACGT ATGT ACTA ATCG TTTT ATTA TTTA"),
        ("T" * 10, "TTTT"),
        ("ACGTCGTAAAA", "TCGT"), # partial match can shadow an exact match
        ]

    nfailed = 0
    for genome, patterns in tests:
        print "genome:", genome
        for pattern in patterns.split():
            print pattern
            a1 = list(WM_approx_Ham1_search(pattern, genome))
            a2 = list(long_check(pattern, genome))
            a3 = list(Paul_McGuire_regex(pattern, genome))
            print a1
            print a2
            print a3
            print a1 == a2, a2 == a3
            nfailed += (a1 != a2 or a2 != a3)
    print "***", nfailed