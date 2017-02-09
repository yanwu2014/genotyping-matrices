### Module with functions for collapsing barcodes within a given edit/hamming distance ###

from editdistance import eval
from collections import defaultdict,Counter
from random import sample
from operator import itemgetter
from time import time
from dedup_umi import ClusterAndReducer

# Input: dictionary with strings as keys, either Counter or int as values; edit distance
# Output: dictionary with keys within edit_dist collapsed, vals summed
def collapseBarcodesExact(barcodeSpacerCounts, edit_dist, barcodeCounts=None, 
                          hamming = False):
    cs = ClusterAndReducer()
    dataType = set

    if not barcodeCounts:
        barcodeCounts = barcodeSpacerCounts
        dataType = int
    else:
        assert set(barcodeCounts.keys()) == set(barcodeSpacerCounts.keys())
    
    if edit_dist == 0 and dataType == int:
        return barcodeSpacerCounts
    elif edit_dist == 0 and dataType != int:
        return barcodeSpacerCounts,barcodeCounts
    
    adj_list = cs.get_adj_list(barcodeCounts.keys(), barcodeCounts, edit_dist, hamming)
    clusters = cs.get_connected_components(barcodeCounts.keys(), adj_list, barcodeCounts)
    finalBarcodes,barcodeCounts = cs.reduce_clusters(clusters, adj_list, barcodeCounts)

    if dataType == int:
        return barcodeCounts
    else:
        newBarcodeSpacerCounts = defaultdict(Counter)
        for parentBarcode,cluster in finalBarcodes.items():
            newCounts = Counter()
            for barcode in cluster:
                newCounts += barcodeSpacerCounts[barcode] 
            newBarcodeSpacerCounts[parentBarcode] = newCounts
    
        return newBarcodeSpacerCounts,barcodeCounts

# Input: iterable, string to query, edit distance 
# Output: item in fuzzySet within edit_dist of value if it exists
#         otherwise return False
def exactExists(myIterable, myStr, edit_dist = 1):
    if myStr in myIterable:
        return myStr
    exists = [(item,eval(item,myStr)) for item in myIterable \
              if eval(item,myStr) <= edit_dist]

    if len(exists) == 0: return False
    elif len(exists) == 1: return exists[0][0]
    else:
        min_dist = min([x[1] for x in exists])
        exists = [x[0] for x in exists if x[1] == min_dist]
        if len(exists) == 1:
            return exists[0]
        else:
            return exists

# Input: dictionary
# Output: key data type, value data type
def checkDataType(myDict):
    k,v = myDict.popitem()
    myDict[k] = v
    return (type(k),type(v))

# Input: vector of ints, minScore
# Output: number of ints in quals < minScore
def belowQual(quals,minScore):
    return sum([x < minScore for x in quals])

# Input: list of tuples (score,item)
# Output: item of with the smallest score
def argmin(tupleList, scorePos=0):
    if scorePos == 0:
        return min(tupleList, key = itemgetter(0))[1]
    if scorePos == 1:
        return min(tupleList, key = itemgetter(1))[0]

# Input: list of tuples (score,item)
# Output: item of with the largest score
def argmax(tupleList, scorePos=0):
    if scorePos == 0:
        return max(tupleList, key = itemgetter(0))[1]
    if scorePos == 1:
        return max(tupleList, key = itemgetter(1))[0]

# Input: string s
# Output: True if s can be converted to int, False otherwise
def isInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


