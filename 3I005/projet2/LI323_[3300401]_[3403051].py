# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 16:26:11 2016

@author: 3300401
"""

import re
import numpy as np
import math
import matplotlib.pyplot as plt
import collections
import heapq


def lireDtrain(filename):
    f = open(filename, "r")
    tabDtrain = list();    
    for line in f:
        if line[0] != '>':
            tmp = re.sub("\s","",line)
            tabDtrain.append(tmp)
    f.close()
    return tabDtrain

def lireTestSeq(filename):
    f = open(filename, "r")
    tabSeq = list()
    tmp = None
    for line in f:
        if line[0] == '>':
            if tmp != None:
                tmp = re.sub("\s","",tmp)
                tabSeq.append(tmp)
            tmp = ""
        else:
            tmp += line
    tmp = re.sub("\s", "", tmp)
    tabSeq.append(tmp)
    f.close()
    return tabSeq
    
def occuDtrain(tabDtrain):
    alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    tabDicoDtrainOccu = list()
    ligne = len(tabDtrain)
    colonne = len(tabDtrain[0])
    for j in range(colonne):
        dico = {cle:0 for cle in alphabet}
        for i in range(ligne):
            dico[tabDtrain[i][j]] += 1
        tabDicoDtrainOccu.append(dico)
    return tabDicoDtrainOccu
    
def poidDtrain(tabDicoDtrainOccu):
    q = 21 # #A  
    tabDicoDtrainPoid = list()
    for i in range(len(tabDicoDtrainOccu)):
        dic = tabDicoDtrainOccu[i]
        nb = np.sum(list(dic.values())) + q
        dico = {cle: (val+1.0)/nb for cle,val in dic.items()}
        tabDicoDtrainPoid.append(dico)
    return tabDicoDtrainPoid
    
def listEntropie(tabDicoDtrainPoid):
    q = 21 # #A
    tabE = list()
    for i in range(len(tabDicoDtrainPoid)):
        entropie = math.log(q, 2)
        entropie += np.sum([poid*math.log(poid,2) for poid in tabDicoDtrainPoid[i].values()])
        tabE.append(entropie)
    return tabE
    
def drainEntropie(tabE):
    plt.plot(tabE)
    plt.xlabel("numero de colonne")
    plt.ylabel("valeur d'entropie relative")
    plt.title("entropie relative en fonction de colonne")
    plt.show()
    
def getMaxPSWM(tabDicoDtrainPoid, i):
    d = list(tabDicoDtrainPoid[i].items())
    d.sort(key = lambda x:x[1], reverse=True)
    return d[0][0]
    
def nBestPosCons(tabDicoDtrainPoid, tabE, n):
    d = [(i,tabE[i]) for i in range(len(tabE))]
    d.sort(key = lambda x:x[1], reverse=True)
    l = list()
    for i in range(n):
        col = d[i][0]
        l.append((col, getMaxPSWM(tabDicoDtrainPoid, col)))
    return l
    
def probInDtrain(tabDicoDtrainPoid, lettre):
    sum = 0
    for i in range(len(tabDicoDtrainPoid)):
        sum += tabDicoDtrainPoid[i][lettre]
    return sum/len(tabDicoDtrainPoid)
    
def modelNull(seq, tabDicoDtrainPoid):
    prob = 1;
    for lettre in seq:
        prob *= probInDtrain(tabDicoDtrainPoid, lettre)
    return prob
    
def logVraisemblance(seq, tabDicoDtrainPoid):
    l = 0
    for i in range(len(seq)):
        l += math.log(tabDicoDtrainPoid[i][seq[i]]/probInDtrain(tabDicoDtrainPoid, seq[i]),2)
    return l
    
def tabLogVarSeq(seq, tabDicoDtrainPoid):
    tablog = list()
    l = len(tabDicoDtrainPoid)
    for i in range(len(seq)-l):
        tablog.append(logVraisemblance(seq[i:i+l], tabDicoDtrainPoid))
    return tablog
    
def drainLogVar(tabLogV):
    plt.plot(tabLogV)
    plt.xlabel("position de depart")
    plt.ylabel("valeur de log vraisemblqnce")
    plt.title("log vraisemblance en fonction du position de depart")
    plt.show()
    
#VmatDtrain = lireDtrain("Dtrain.txt")
#VmatTestSeq = lireTestSeq("test_seq.txt")
#VtabDicoDtrainOccu = occuDtrain(VmatDtrain)
#VtabDicoDtrainPoid = poidDtrain(VtabDicoDtrainOccu)
#VtabE = listEntropie(VtabDicoDtrainPoid)
#VtreeBestPosCons = nBestPosCons(VtabDicoDtrainPoid, VtabE, 3)
#VtabLogV = tabLogVarSeq(VmatTestSeq[0], VtabDicoDtrainPoid)

#print(VmatDtrain)
#print(VmatTestSeq)
#print(VtabDicoDtrainOccu)
#print(VtabDicoDtrainPoid)
#print(VtabE)
#drainEntropie(VtabE)
#print(VtreeBestPosCons)
#print(VtabDicoDtrainPoid[0]['-']) #should be equal to 0.31
#print(VtabE[0]) #should be equal to 1.85
#print(VtabLogV)
#drainLogVar(VtabLogV)

def mat_co_occu(matDtrain):
    M = len(matDtrain)
    L = len(matDtrain[0])
    mat = np.zeros((L,L), dtype=collections.Counter)
    for i in range(L):
        for j in range(i,L):
            l = [matDtrain[c][i]+matDtrain[c][j] for c in range(M)]
            mat[i][j] = collections.Counter(l)
    return mat
   
def poid_co_occu(i,j,a,b,mat,M,q):
    return (mat[i][j][a+b]+1.0/q)/(M+q)
    
def inf_Mutuelle(matDtrain, tabDicoDtrainPoid):
    mat = mat_co_occu(matDtrain)
    alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    q = len(alphabet)
    M = len(matDtrain)
    L = len(matDtrain[0])
    mutuelle = np.zeros((L,L))
    for i in range(L):
        for j in range(i+1,L):
            sum = 0
            for a in alphabet:
                for b in alphabet:
                    poid_co_o = poid_co_occu(i,j,a,b,mat, M, q)
                    sum += poid_co_o*math.log(poid_co_o/(tabDicoDtrainPoid[i][a]*tabDicoDtrainPoid[j][b]),2)
            mutuelle[i][j] = sum
    return mutuelle
    
def nBestMutVal(matMut, n):
    max = list()
    L = len(matMut)
    count = 0
    for i in range(L):
        for j in range(i, L):
            if count < n:
                heapq.heappush(max, (matMut[i][j], i, j))
                count += 1
            else:
                if matMut[i][j] > max[0][0]:
                    heapq.heapreplace(max, (matMut[i][j], i, j))
    return sorted(max)
    
def lireDistanceFile(filename, L):
    mat = np.zeros((L,L))
    f = open(filename, "r")
    for line in f:
        line = re.sub("\n", "", line)
        line = re.split("\s", line)
        if len(line) == 3:
            mat[int(line[0])][int(line[1])] = float(line[2])
    return mat
            
def prop(matDistance, matMuttuel, nb, distance):
    count = 0
    tabBest = nBestMutVal(matMuttuel, nb)
    for n in range(nb):
        i = tabBest[n][1]
        j = tabBest[n][2]
        if matDistance[i][j] > 0 and matDistance[i][j] < distance:
            #print(str(i)+" "+str(j)+" "+str(matDistance[i][j]))
            count += 1
    return count/float(nb)
    
def traceFonc(matDistance, matMuttuel, nb, distance):
    l = [0]
    for i in range(1,nb+1):
        l.append(prop(matDistance, matMuttuel, i, distance))
    plt.plot(l)
    plt.axis([0,nb+1,0,1.1])
    plt.xlabel("nombre de plus grande valeur de infrmation mutuelle considerees")
    plt.ylabel("proportion des paires qu ont une distance inferieur a "+str(distance))
    plt.title("proportion des paires de distance < "+str(distance)+" en fonction de nb considerees")
    plt.show()
           
#VmatMutuelle = inf_Mutuelle(VmatDtrain, VtabDicoDtrainPoid)  
#VmatDistance = lireDistanceFile("distances.txt", 48)
#VtabnBestMut = nBestMutVal(VmatMutuelle, 50)
#VtabFaction = prop(VmatDistance, VmatMutuelle, 1, 8)
      
#print(VmatMutuelle[0][1])
#print(VtabnBestMut)
#print(VmatDistance)
#print(VtabFaction)
#traceFonc(VmatDistance, VmatMutuelle, 50, 8)

def outputEntroRel(filename, tabE, tabDicoDtrainPoid):
    f = open(filename, "w")
    for i in range(len(tabE)):
        f.write(str(i)+"\t"+str(tabE[i])+"\t"+str(getMaxPSWM(tabDicoDtrainPoid,i))+"\n")
    f.close()
    
def outputLogV(filename, tabLogV):
    f = open(filename, "w")
    for i in range(len(tabLogV)):
        f.write(str(i)+"\t"+str(tabLogV[i])+"\n")
    f.close()
    
def outputBestMut(filename, tabFaction):
    tabFaction = sorted(tabFaction, reverse=True)
    f = open(filename, "w")
    for i in range(len(tabFaction)):
        f.write(str(tabFaction[i][1])+","+str(tabFaction[i][2])+"\t"+str(tabFaction[i][0])+"\n")
    f.close()
    
if __name__ == "__main__":
    VmatDtrain = lireDtrain("Dtrain.txt")
    VmatTestSeq = lireTestSeq("test_seq.txt")
    VtabDicoDtrainOccu = occuDtrain(VmatDtrain)
    VtabDicoDtrainPoid = poidDtrain(VtabDicoDtrainOccu)
    VtabE = listEntropie(VtabDicoDtrainPoid)
    VtreeBestPosCons = nBestPosCons(VtabDicoDtrainPoid, VtabE, 3)
    VtabLogV = tabLogVarSeq(VmatTestSeq[0], VtabDicoDtrainPoid)
    
    #print(VmatDtrain)
    #print(VmatTestSeq)
    #print(VtabDicoDtrainOccu)
    #print(VtabDicoDtrainPoid)
    #print(VtabE)
    #drainEntropie(VtabE)
    #print(VtreeBestPosCons)
    #print(VtabDicoDtrainPoid[0]['-']) #should be equal to 0.31
    #print(VtabE[0]) #should be equal to 1.85
    #print(VtabLogV)
    #drainLogVar(VtabLogV)
    
    
    VmatMutuelle = inf_Mutuelle(VmatDtrain, VtabDicoDtrainPoid)  
    VmatDistance = lireDistanceFile("distances.txt", 48)
    VtabnBestMut = nBestMutVal(VmatMutuelle, 50)
    VtabFaction = prop(VmatDistance, VmatMutuelle, 1, 8)
          
    #print(VmatMutuelle[0][1])
    #print(VtabnBestMut)
    #print(VmatDistance)
    #print(VtabFaction)
    #traceFonc(VmatDistance, VmatMutuelle, 50, 8)
    
    #les 3 fichiers demender
    outputEntroRel("entropie_relative.txt",VtabE,VtabDicoDtrainPoid)
    outputLogV("log_vraisemblance.txt", VtabLogV)
    outputBestMut("paires_de_positions.txt", VtabnBestMut)