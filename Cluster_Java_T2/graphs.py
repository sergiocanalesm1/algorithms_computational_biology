from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt

def graph_covid():
    mypath = "ProcsData"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for file in onlyfiles:
        if "Overlap" in file:
            continue
        if "COVID" in file:
            f = open(join(mypath,file))
            ksize = file.split(".")[-2].split("_")[1]
            print(ksize)
            ar = list(map(int, f.read().rstrip().split()))
            x = [ar[i] for i in range(len(ar)) if i % 2 == 0]
            y = [ar[i] for i in range(len(ar)) if i % 2 != 0]
            #plt.figure()
            if file.startswith("COVID"):
                profundidad = int(file.split("X")[0].split("fasta")[1])


            """  
            plt.figure(file)
            plt.title(" Kmer abundance for {}".format(file))
            plt.ylabel("Abundance")
            plt.xlabel("# of kmers")
            plt.bar(x,y)
            plt.savefig('graphs/{}.png'.format(file), dpi=300, bbox_inches='tight')
            
            """
            f.close()



def graph_test():
    colors = {
        5:"red",
        10:"blue",
        15:"yellow",
        20:"brown",
        25:"green",
        50:"black",
        75:"pink"
    }
    mypath = "ProcsData"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2,figsize=(10, 8))
    sub = {"10":(ax1,50),
            "20":(ax2,60),
           "50": (ax3,80),
           "100": (ax4,150)
           }
    for file in onlyfiles:
        if "kmer" in file and "COVID" not in file:
            f = open(join(mypath, file))
            ksize = int(file.split(".")[-2].split("_")[1])
            profundidad = file.split("x")[0].split("test")[1]
            ar = list(map(int, f.read().rstrip().split()))
            x = [ar[i] for i in range(len(ar)) if i % 2 == 0]
            y = [ar[i] for i in range(len(ar)) if i % 2 != 0]
            sub[profundidad][0].bar(x,y,color=colors[ksize],label=ksize)
            f.close()
    for key, value in sub.items():
        value[0].set_xlim([0,value[1]])
        value[0].set_ylabel("# of kmers")
        value[0].set_xlabel("Abundance")
        value[0].title.set_text("profundidad: {}x".format(key))

    handles, labels =sub["10"][0].get_legend_handles_labels()
    fig.legend(handles, labels,title="Kmer size",loc="center right")
    #plt.xlabel("# of kmers")
    #plt.ylabel("Abundance")
    plt.subplots_adjust(hspace=0.8,wspace=0.5)
    fig.savefig('graphs/{}.png'.format(file), dpi=300)

def mean():
    mypath = "ProcsData"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for file in onlyfiles:
        if "test" in file and "kmer" in file:
            f = open(join(mypath, file))
            ar = list(map(int, f.read().rstrip().split()))
            x = [ar[i]* i for i in range(len(ar)) if i % 2 == 0]
            print("file: {} mean kmer abundance: {}".format(file,sum(x)))

def graph_5():
    mypath = "ErrorAnalisis"
    onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(8, 4))
    for file in onlyfiles:
        if file.endswith(".dat"):
            f = open(join(mypath,file))
            error = int(file.split(".")[-2].split("Error")[1])
            ar = list(map(int, f.read().rstrip().split()))
            x = [ar[i] for i in range(len(ar)) if i % 2 == 0]
            y = [ar[i] for i in range(len(ar)) if i % 2 != 0]
            ax1.bar(x,y,label="{}%".format(error),alpha=0.4)
            if error == 0:
                ax2.set_ylabel("# of kmers")
                ax2.set_xlabel("Abundances")
                ax2.bar(x, y,color="#FE8108", label="{}%".format(error),alpha=0.4)
                ax2.legend(title="Error")
            f.close()
    ax1.set_yscale('log')
    ax1.set_ylabel("# of kmers")
    ax1.set_xlabel("Abundances")
    fig.suptitle("Kmer distribution on sequence length 50 with mean 20x")
    ax1.legend(title="Error")
    plt.subplots_adjust( wspace=0.5)
    #fig.savefig('{}/graph_analisis.png'.format(mypath), dpi=300)
    plt.show()

mean()
