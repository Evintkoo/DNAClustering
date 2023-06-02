import numpy as np
import math
import random
import sys
import pandas as pd
import time
import csv
import multiprocessing as mp
from multiprocessing import Pool

'''generateDNA: generate a random DNA strain with length "length"'''
# returns DNA strain
def generateDNA(length: int) -> str:
    # the chromosome is either "a", "g", "t", "c"
    choice = "agtc"
    
    # pick random of chromosome and arrange into a dna strain
    s = "".join(random.choice(choice) for x in range(length))
    return s

'''CLASS OF DNA MUTATION: DELETION, DUPLICATION, INVERSION, INSERTION, TRANSLOCATION, POINT MUTATION, AND FRAMESHIFT '''

class mutate:
    
    '''deletion: delete random number of random chromosome in DNA'''
    # mutate.deletion() giving an input of string DNA strain and remove several DNA strain randomly
    # returns deleted DNA strain 
    def deletion(DNAStrain: str) -> str:
        # select a number of deleted genoms
        j = random.randint(0,len(DNAStrain)-1)
        s =  DNAStrain
        
        # deleting random character in DNAStrain with total of j deletion
        for delete in range(0,j):
            # select a random index of string in DNA strain
            i = random.randint(0,len(s)-1)
            
            # delete an character in the index
            s =  s[:i] + s[i+1:]
    
        return s
    
    '''duplication: duplicate 1 random DNA segment'''
    # mutate.duplication() giving an input of string DNA strain and duplicate a segment of the DNA strain with random length and random place
    # return a DNA strain of duplicated DNA segment
    def duplication(DNAStrain: str) -> str:
        # select a random left index of DNA segment
        left = random.randint(0,len(DNAStrain)-1)
        
        # select a random right index of DNA segment
        right = random.randint(left, len(DNAStrain)-1)
        
        # duplicate an character in the index
        s =  DNAStrain[:left] + DNAStrain[left:right] + DNAStrain[left:right] + DNAStrain[right:]
    
        return s
    
    '''inversion: insert a random DNA segment to the strain'''
    # mutate.inversion() giving an input of string DNA strain and duplicate 1 of the DNA strain randomly
    # returns inverted DNA strain
    def inversion(DNAStrain: str) -> str:
        # select a random left index of DNA segment
        left = random.randint(0,len(DNAStrain)-1)
        
        # select a random right index of DNA segment
        right = random.randint(left, len(DNAStrain)-1)
        
        # take a part of strain in DNA and reverse it
        s = DNAStrain[left:right]
        s = s[::-1]
        
        # delete an character in the index
        s =  DNAStrain[:left] + s + DNAStrain[right+1:]
    
        return s
    
    '''insertion: insert a random DNA segment'''
    # mutate.insertion() giving an input of string DNA strain and insert a DNA segment into DNA strain randomly
    # returns inserted DNA strain
    def insertion(DNAStrain: str) -> str:
        # select a random length of string of inserted DNA segment
        length = random.randint(0,len(DNAStrain)-1)
        
        # generate random DNA segment
        s = generateDNA(length)
        
        # select a random index to insert the DNA segment
        i = random.randint(0,len(DNAStrain)-1)
        
        # delete an character in the index
        s =  DNAStrain[:i] + s + DNAStrain[i+1:]
    
        return s
    
    '''translocation: move a random DNA segment'''
    # mutate.inversion() giving an input of string DNA strain and duplicate 1 of the DNA strain randomly
    # returns translocated DNA
    def translocation(DNAStrain: str) -> str:
        # select a random index of string in DNA strain
        i = random.randint(0,len(DNAStrain)-1)
        
        # select a random length of string of inserted DNA segment
        length = random.randint(0,len(DNAStrain)-1)
        
        # generate random DNA segment
        s = generateDNA(length)
        
        # delete an character in the index
        s =  DNAStrain[:i] + s
    
        return s
    
    '''Point mutation: change a single nucleotide in DNA'''
    # mutate.pointMutation giving an input of string of DNA strain and change a char inside it
    # returns modified DNA
    def pointMutation(DNAStrain: str) -> str:
        # select a random index of string in DNA strain
        i = random.randint(0,len(DNAStrain)-1)
        
        # generate a random nucleotide
        s = generateDNA(1)
        
        # delete an character in the index
        s =  DNAStrain[:i] + s + DNAStrain[i+1:]
    
        return s
    
    '''Frameshift mutation: insert a single nucleotide in DNA'''
    # mutate.frameshift giving an input of string of DNA strain and insert a random nucleotide in it
    # return modified DNA strain
    def frameshift(DNAStrain: str) -> str:
        # select a random index of string in DNA strain
        i = random.randint(0,len(DNAStrain)-1)
        
        # generate a random nucleotide
        s = generateDNA(1)
        
        # delete an character in the index
        s =  DNAStrain[:i] + s + DNAStrain[i:]
    
        return s

# strainDifference inputing 2 string and count the difference
# return a float [0,1] that represent the similarity between 2 string (1: similar) (0: totally different)
def strainSimilarity (str1: str, str2: str) -> float:
    str1 = str1 + (' ' * (len(str2) - len(str1)))
    str2 = str2 + (' ' * (len(str1) - len(str2)))
    return sum(1 if i == j else 0
               for i, j in zip(str1, str2)) / float(len(str1))

# define if the mutation is a cancer or not based on the difference of the string
# return true if the mutation is a cancer, else return false
def isCancer(strain1: str,strain2: str) -> bool:
    dif = strainSimilarity(strain1,strain2) 
    
    if dif < cancerRate:
        return True
    else :
        return False


'''Class of the agent which save the information of its ID: the ID of the agent (int), 
                DNAstrain: unique DNA string of agent (string), and cancerstatus: a boolean that represent the status of agent if it has cancer or not (boolean)'''
# agent is contains:
# ID: unique id of each agent, 
# generationNumber: represents the generation that agent have, 
# DNAstrain: the DNA of agent, 
# status: represent if its a cancer or not
class agent:
    def __init__ (self, id: str, generationNumber: int, DNAstrain: str, status: bool):
        self.id = id
        self.generationNumber= generationNumber
        self.DNAstrain = DNAstrain
        self.status = status


'''Class of a tree data structure that contain the agents'''
class TreeNode:
    # for each TreeNode, contains a data of agent, list of children, and the parent
    def __init__(self, data: agent):
        self.data = data
        self.children = []
        self.parent = None
    
    # adding a child on the TreeNode
    def addChild (self, child):
        child.parent = self
        self.children.append(child)
    
    # return the level/generation of the data
    def getLevel(self):
        level = 0
        p = self.parent
        while p:
            level += 1
            p = p.parent
        
        return level
    
    # print the tree structure
    def printTree(self):
        print("/t" * self.getLevel(), [self.data.id, self.data.generationNumber, self.data.DNAstrain, self.data.status])
        if self.children:
            for child in self.children:
                child.printTree()
    
    # return list of agent that has cancer
    def cancerNode(self, cancerList: list) -> list:
        if self.data.status == True:
            cancerList.append(self.data)
            # print([self.data.id, self.data.status])
        if self.children:
            for child in self.children:
                child.cancerNode(cancerList)
        return cancerList
    
    # save the tree in a file
    # make a data, save each node as:
    # id, generationNumber, DNAstrain, status, parent.id
    def dataTree(self,data: list) -> list:
        if self.parent: 
            data.append([self.data.id, self.data.generationNumber, self.data.DNAstrain, self.data.status, self.parent.data.id])
        if self.children:
            for childs in self.children:
                childs.dataTree(data)
        return data
        

'''isMutated() is decide if the mutation will happen or not'''
# return true if it will happen, and return false if not
def isMutated() -> bool:
    x = random.uniform(0,1)
    if x < mutationRate:
        return True
    else :
        return False
        


'''pickMutation is a function that will randomly choose the type of mutation and implement the picked mutation into the DNA strain'''
# return the mutated DNA
def pickMutation(DNAstrain: str) -> str:
    rand = random.choice(range(0,7)) # weighting
    match rand:
        case 0:
            r = mutate.deletion(DNAstrain)
            pass
        case 1:
            r = mutate.duplication(DNAstrain)
            pass
        case 2:
            r = mutate.inversion(DNAstrain)
            pass
        case 3:
            r = mutate.insertion(DNAstrain)
            pass
        case 4:
            r = mutate.translocation(DNAstrain)
            pass
        case 5:
            r = mutate.pointMutation(DNAstrain)
            pass
        case 6:
            r = mutate.frameshift(DNAstrain)
            pass
    return r

'''idMaker is a function that make a unique ID for each agent'''
# returns string of number that represent the unique ID
def idMaker(idList: list) -> str:
    num = "0123456789abcdefghijklmnopqrstuvxyz"
    idLength = 20
    id = "".join(random.choices(num, k=idLength))
    while id in idList:
        id = "".join(random.choices(num, k=idLength))
    listID.append(id)
    # print(id)
    idList = idList.sort()
    return id            
    
'''mutation is a function that generate the mutation of the agent'''
# return agent with new ID, DNAstrain
def mutation(parentAgent: agent, genNumber: int) -> agent:
    if isMutated():
        childDNAMutation = pickMutation(parentAgent.DNAstrain)
        return agent(id= idMaker(listID),
                    generationNumber=genNumber, 
                    DNAstrain=childDNAMutation, 
                    status=isCancer(parentAgent.DNAstrain,childDNAMutation))
    else: 
        return agent(id = idMaker(listID),
                    generationNumber=genNumber, 
                    DNAstrain = parentAgent.DNAstrain, 
                    status = False)
    

'''buildGenerationTree is a function build a tree data structure of the generation'''
# returns tree of genetical mutation
def buildGenerationTree(parentAgent: agent, generation: int, totalAgent: int) -> TreeNode:
    totalAgent = totalAgent+1
    # print([parentAgent.id])
    parrent = TreeNode(parentAgent)
    # print(parentAgent.id)
    
    
    if generation > 0 and parentAgent.status == False:
        # make a mutation of next generation
        for i in range(Numchild):
            child = mutation(parentAgent, totalAgent)

            parrent.addChild(buildGenerationTree(parentAgent=child, 
                            generation= generation-1,
                            totalAgent= totalAgent))
    return parrent


'''mutationPath is a function that finds the acentors path of the agent based on its ID'''
# returns a list of ID from parentDNA to the agent
def mutationPath(env: TreeNode, id: str, path: list) -> list:
    #print(path)
    if env.data.id == id:
        path.append(env.data)
        return path
    elif env.children: 
        path.append(env.data)
        for childs in env.children:
            ret = mutationPath(childs, id, path)
            if ret:
                return ret
        path.remove(env.data)


'''findAgent is a function that finds the acentors path of the agent based on its ID'''
# returns a agent object of targeted ID
def findAgent(env: TreeNode, id: str) -> agent:
    if env.data.id == id:
        return env.data
    else :
        for childs in env.children:
            ret = findAgent(childs, id)
            if ret:
                return ret

# saveTree is save the tree into CSV data 
def saveTree(env: TreeNode, path: str) -> None:
    table = environment.dataTree([])
    columns = ["Agent ID", "Generation Number", "DNA Strain", "Cancer Status", "Parent ID"]
    df = pd.DataFrame(table, 
                      columns=columns)
    df.to_csv(path, index=False)

# check the similarity based on its nucleotids
def DNASimilarity (DNAStrain1: str, DNAStrain2: str) -> float:
    ret = 0
    for i in range (min(len(DNAStrain1), len(DNAStrain2))):
        if DNAStrain2[i] == DNAStrain1[i]:
            ret += 1
    return round(ret/min(len(DNAStrain1), len(DNAStrain2))*100,2)

def build3DdnaMatrix(listDNA: list[str])  :
    matrix = []
    for i in range(len(listDNA)):
        matrixX = []
        for j in range(len(listDNA)):
            x = DNASimilarity(listDNA[i].DNAstrain, listDNA[j].DNAstrain)
            matrixX.append(x)
        matrix.append(matrixX)
    return matrix
            

if __name__ == "__main__" :
    
    '''PREDEFINED VARIABLE'''
    start_time = time.time()
    #similarity of DNA strain to be define as cancer
    cancerRate = strainSimilarity("gaaacctgtttgttggacatactggatacagctggacaagaagagtacagtgccatgagagaccaatacatgaggacaggcgaaggcttcctctgtgtatttgccatcaataat",
                                  "gaaacctgtttgttggatactttgcttttgctcactcgaaggcttcctctgtgtatttgccatcaataa")      
    # maximum number of parent to have the child on the next gen
    Numchild = 2          
    # 1 in 10.000 generation is mutating   
    mutationRate = 0.06
    #the original DNA
    parentDNA = []
    print("--- Load the DNA datasets ---")
    for i in range(1,26):
        with open('./Parrent Database/{}.txt'.format(i), 'r') as f:
            for line in f.readlines():
                l = line.replace(" ", "")
                parentDNA.append(l)
    print("Total DNA parent:", len(parentDNA))
    
    #total repetition per DNA parent
    totalRepetition = 50
    
    # total generation generated
    totalGeneration = 14
    
    #empty list of ID
    listID = []
    
    #setting recusion limit to the device
    sys.setrecursionlimit(100000)
    
    table = []
    #start the simulation
    for DNAid in range(len(parentDNA)):
        print("Parent DNA", DNAid, "---------------------------------------------------")
        for reps in range(totalRepetition):
            print("Repetition", reps+1, "----------------------------------")
            # generate a parent agent
            parental = agent(id=idMaker(listID),
                            generationNumber=0, 
                            DNAstrain=parentDNA[DNAid], 
                            status=False)
            
            #build the tree data structure
            print("---- generating the environment --")
            environment = buildGenerationTree(parental, totalGeneration, 0)
            print("---- done building environment ---")
            print("total agent generated:", len(listID))
            #save the environment data
            #path = "./Datas/AgentsData{}_{}.csv".format(DNAid, reps)
            #print("saving tree to", path)
            #saveTree(env=environment, 
            #        path= path)
            
            #make a list of cancer agent from the tree
            cancerAgent = environment.cancerNode(cancerList=[])
            print("total cancer agents:", len(cancerAgent))
            
            #make a list of path that leads to cancer
            cancerPathAgentList = []
            for agents in cancerAgent:
                x = mutationPath(env= environment,
                                    id= agents.id,
                                    path=[])
                cancerPathAgentList.append(x)
                #print(x)
            
            #make a list of agent that leads to cancer mutation
            print("--- generating path of cancer-leading agent ---")
            listAgentLeadToCancer = []
            for agentPaths in cancerPathAgentList:
                similarityList = []
                for i in range(len(agentPaths)-2):
                    similarityList.append(strainSimilarity(agentPaths[i].DNAstrain,agentPaths[i+1].DNAstrain))
                # print(similarityList)
                mins = 1
                if similarityList:
                    mins = min(similarityList)
                if mins < 1:
                    minIndex = similarityList.index(mins)
                    listAgentLeadToCancer.append(agentPaths[minIndex+1])
            print("total agent leading to cancer:", len(listAgentLeadToCancer))
            for agents in listAgentLeadToCancer:
                table.append([agents.id, agents.generationNumber, agents.DNAstrain, agents.status])
            
            #building the 3D matrix clusterring
            #print("--- building the matrix clustering --")
            #matrix3D = build3DdnaMatrix(listAgentLeadToCancer)
            
            #save the cancer-leading agent
            #columns = [agents.id for agents in listAgentLeadToCancer]
            #df = pd.DataFrame(matrix3D, 
            #                columns=columns)
            #df.to_csv("./Datas/similarityMatrix{}_{}.csv".format(i, j), index=False)
            listID = []
    
    columns = ["Agent ID", "Generation Number", "DNA Strain", "Cancer Status"]
    df = pd.DataFrame(table, 
                    columns=columns)
    df.to_csv("./Datas/agentLeadToCancer.csv", index=False)    
    print("--- end of the program ---")
    print("Program executed in %s seconds " % (time.time() - start_time))
    