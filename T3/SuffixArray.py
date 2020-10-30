import sys
import numpy as np
s=sys.argv[1]
#s="Banana$"
from collections import Counter
SuffixArray=([t[1] for t in sorted((s[i:],i) for i in range(len(s)))])
print("Suffix array")
print(SuffixArray)
print("Burrows-Wheeler matrix")
for i in SuffixArray:
    print(s[i:]+s[:i])
print("Burros-Wheeler transform")
BWT=[]
for i in SuffixArray:
    print(s[i-1],end="")
    BWT.append(s[i-1])
print()
Counts=dict(Counter(BWT))
Tally=np.zeros((len(s)+1,len(Counts)+1))
print(BWT)
