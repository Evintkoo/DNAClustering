listDF = list()
    for dicts in kmerList:
        dictToList = list()
        for keys in X:
            try:
                dictToList.append(dicts[keys])
            except:
                dictToList.append(0)
        listDF.append(dictToList)
    y = model.predict(listDF)