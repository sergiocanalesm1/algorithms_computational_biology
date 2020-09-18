from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt

def graph():
    mypath = "ProcsData"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for file in onlyfiles:
        if "Overlap" in file:
            continue
        f = open(join(mypath,file))
        #ksize = file.split("_")[-1]
        ar = list(map(int, f.read().rstrip().split()))
        x = [ar[i] for i in range(len(ar)) if i % 2 == 0]
        y = [ar[i] for i in range(len(ar)) if i % 2 != 0]
        plt.figure(file)
        plt.title(" Kmer abundance for {}".format(file))
        plt.xlabel("Abundance")
        plt.ylabel("# of kmers")
        plt.bar(x,y)
        plt.savefig('graphs/{}.png'.format(file), dpi=300, bbox_inches='tight')
        f.close()

graph()


