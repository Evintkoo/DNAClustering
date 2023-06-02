import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import math
import sklearn
import pickle
from numpy import unique
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
import matplotlib.pyplot as plt

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

if __name__ == "__main__":
    print("Load the model")
    with open('C:/Users/Evint/Documents/College Classes/VANT 148 VGM/Source code/models/DBScan_model.pkl', 'rb') as f:
        model = pickle.load(f)
    path = "C:/Users/Evint/Documents/College Classes/VANT 148 VGM/Source code/Datas/agentLeadToCancer.csv"
    df = pd.read_csv(path)
    listDNA = list(df["DNA Strain"])
    
    print("load the clusters data")
    path = "C:/Users/Evint/Documents/College Classes/VANT 148 VGM/Source code/Datas/clusters.csv"
    df_cluster = pd.read_csv(path)
    kmerList = list()
    for strains in listDNA:
        kmerList.append(count_kmers(strains,int(len(df_cluster.columns[0]))))
    len(kmerList)
    X = df_cluster.columns
    print("cleaning data")
    listDF = list()
    for dicts in kmerList:
        dictToList = list()
        for keys in X:
            try:
                dictToList.append(dicts[keys])
            except:
                dictToList.append(0)
        listDF.append(dictToList)
    y = model.fit_predict(listDF)
    
    print("clasifying data")
    dtree = DecisionTreeClassifier(max_depth=5)
    dtree = dtree.fit(listDF, y)
    y = [str(x) for x in y]
    
    print("plot data")
    plt.figure(figsize=(650,100))
    tree.plot_tree(dtree, 
                   feature_names=X, 
                   class_names=unique(y),
                   fontsize=5, 
                   max_depth=5, 
                   impurity=False,
                   precision=2,
                   filled=True)
    
    print("show plot")
    plt.show()
    
    data = pd.DataFrame(listDF, columns=df_cluster.columns)
    data["cluster"] = y
    data.to_csv("Datas/clusteredAgentLeadToCancer.csv")
    