# Libraries
import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import math
from sklearn import datasets
import pickle
import time
from sklearn_som.som import SOM
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree

def count_kmers(sequence, k_size):
    data = {}
    size = len(sequence)
    for i in range(size - k_size + 1):
        kmer = sequence[i: i + k_size]
        try:
            data[kmer] += 1
        except KeyError:
            data[kmer] = 1
    return data

def Merge(dict1, dict2):
    return(dict2.update(dict1))

def DNAClustering(listDNA, m, n):
    print("Encode the DNA list")
    averagelength = sum([len(item) for item in listDNA])/len(listDNA)
    int(math.log(averagelength,4))
    kmerList = list()
    for strains in listDNA:
        kmerList.append(count_kmers(strains,int(math.log(averagelength,4))))
        

    X = {}
    print("encode the parameter of the list")
    for dicts in kmerList:
        X.update(dicts)
    X = list(X.keys())
    params = X
    listDF = list()
    for dicts in kmerList:
        dictToList = list()
        for keys in X:
            try:
                dictToList.append(dicts[keys])
            except:
                dictToList.append(0)
        listDF.append(dictToList)
    
    print("start clustering")
    
    clustering = SOM(m=m, n=n).fit(listDF)
    return [params]

if __name__ == "__main__" :
    sys.setrecursionlimit(100000)
    start_time = time.time()
    print("extracting data ...")
    path = "C:/Users/Evint/Documents/College Classes/VANT 148 VGM/Source code/Datas/agentLeadToCancer.csv"
    print("data from {}".format(path))
    df = pd.read_csv(path)
    df = df.drop_duplicates()
    listDNA = df["DNA Strain"].values
    print("total DNA strains:", len(listDNA))
    a = DNAClustering(listDNA, 10, 10)
    #clusters = pd.DataFrame(a[1], columns= a[0])
    #clusters.to_csv("C:/Users/Evint/Documents/College Classes/VANT 148 VGM/Source code/Datas/clusters.csv", index=False)
    print(a)
    print("Program executed in %s seconds " % (time.time() - start_time))