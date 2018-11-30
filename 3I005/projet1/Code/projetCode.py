# -*- coding: utf-8 -*-
"""
Created on Sun Feb  7 14:20:50 2016

@author: Johnny
"""
import numpy as np
import math
import matplotlib.pyplot as plt
import unicodedata
import re
import codecs
import collections
import heapq
import random

p = [1,2,3,4,5,6,7,8,9]
q = [9,8,7,6,5,4,3,2,1]

def normaliserDict(dict):
    cle = list(dict.keys())
    sum = np.sum(list(dict.values()))
    for c in cle:
        dict[c] /= sum
    return dict

def openfile(filename):
    return codecs.open(filename, encoding="utf-8").read()

def entropie(l):
    return -np.sum([0 if x==0 else x*math.log(x,2) for x in l])

def divergence(p, q):
    return np.sum([float("inf") if q[i] == 0 else p[i]*math.log(p[i]/q[i], 2) for i in range(len(p))])

def divergenceDict(p, q):
    p_keys = list(p.keys())
    q_keys = list(q.keys())
    return np.sum([(p[key]*math.log(p[key]/q[key], 2) if key in q_keys else float("inf")) for key in p_keys])
    
def graph():
    x = np.arange(0, 1.01, 0.01)
    y = [entropie([p, 1-p]) for p in x]
    plt.plot(x, y)
    plt.xlabel("valeur p")
    plt.ylabel("entropie")
    plt.title("entropie du loi binomiale de parametre p")

def preprocess(text):
    data = unicodedata.normalize('NFKD', text).encode('ASCII', 'ignore').decode("utf-8").lower()
    data = re.sub("\s", "", re.sub("\s(?=\s)","",data))
    #print(data)
    return data

def count_ngrams(text, n):
    #print(text)
    data = preprocess(text)
    #print(data)
    cou = collections.Counter()    
    for i in range(len(data)-n+1):
        if re.search("[^a-z]",data[i:i+n]) == None:
            cou[data[i:i+n]] += 1
    dic = dict(cou)
    dic = collections.OrderedDict(sorted(dic.items()))
    return dic
    
def entropieDict(dico):
    entrop = entropie(list(dico.values()))
    return entrop
    
def entropieText(text, n):
    dic = normaliserDict(count_ngrams(text, n))
    return entropieDict(dic)
    
def divergenceText(text1, text2, n):
    dico1 = normaliserDict(count_ngrams(text1, n))
    dico2 = normaliserDict(count_ngrams(text2, n))
    return divergenceDict(dico1, dico2)
  
def histograme(text, n, titre="histograme de frequence des lettre"):
    dic = normaliserDict(count_ngrams(text, n))
    index = dic.keys()
    freq = list(dic.values())
    plt.bar(np.arange(len(index)), freq, 1)
    ax = plt.axes()
    ax.set_xticks(np.arange(len(index)) + (1 / 2))
    ax.set_xticklabels(index)
    plt.title(titre+" N="+str(n))
    
def histogrameFile(filename, n):
    histograme(openfile(filename), n, filename)  
    
def deterLangByDict(dico, dicoDeDico):
    val = float('inf')
    lan = ''
    for langue in dicoDeDico.keys():
        div = divergenceDict(dico, dicoDeDico[langue])
        print("val relative de "+langue+" est "+str(div))
        if  div < val:
            val = div
            lan = langue
    return val, lan
#cherchant la valuer d entropie relative minimum
def deterLangEntropie(filename, listfilenameDictio, n):
    dico = normaliserDict(count_ngrams(openfile(filename), n))
    dicoDeDico = dict();
    for name in listfilenameDictio:
        dicoDeDico[name] = normaliserDict(count_ngrams(openfile(name), n))
    divVal, lang = deterLangByDict(dico, dicoDeDico)
    print(filename+' appartient au '+lang+' avec comme entropie relative '+str(divVal))

#dicoDepSup nbDep+1 lettre donc nbDep+1_grams
#dicoDepInf nbDep lettre donc nbDep_grams
def probWL(tabLettre, dicoDepSup, dicoDepInf, nbDep):
    p = 0
    keySup = dicoDepSup.keys()
    #si la proba les 0 alors on fait -100 
    #qui correspond a un proba de 7.88e-31 donc pratiquement 0
    #le but de ne pas utiliser -Infinie est de pouvoir encore comparer
    #meme si dans le dico de la langue correspondant apparait des proba de 0
    if nbDep == 0:#dans le cas ou ya pas de dependance
        for i in range(0, len(tabLettre)):
            p += math.log(dicoDepSup[tabLettre[i]]) if tabLettre[i] in keySup else -100
        return p
    
    keyInf = dicoDepInf.keys()
    #proba de x0=w0...xnbDep=wnbDep
    motinf = ''.join(tabLettre[0:nbDep])
    
    p += math.log(dicoDepInf[motinf], 2) if motinf in keyInf else -100
    for i in range(nbDep, len(tabLettre)):
        motsup = ''
        motinf = ''
        for l in range(i-nbDep, i):
            motsup += tabLettre[l]
            motinf += tabLettre[l]
        motsup += tabLettre[i]
        p += (math.log(dicoDepSup[motsup]/dicoDepInf[motinf] ,2)) if (motinf in keyInf and motsup in keySup) else -100
    return p
    
def createTabLettre(text, n):
    tabLettre = list()
    data = preprocess(text)
    for i in range(len(data)-n+1):
        if re.search("[^a-z]",data[i:i+n]) == None:
            tabLettre.append(data[i:i+n])
    return tabLettre

#nbDep nombre de dependance
#n la taille d une lettre i.e nb de lettre pour constituer un lettre de la formule
def deterLangBaye(filename, listfilenameDictio, nbDep, n):
    max = float('-inf')
    lang = ''
    tabLettre = createTabLettre(openfile(filename), 1)
    for name in listfilenameDictio:
        dicoDepSup = normaliserDict(count_ngrams(openfile(name), nbDep*n+1))
        dicoDepInf = None
        if nbDep > 0:
            dicoDepInf = normaliserDict(count_ngrams(openfile(name), nbDep*n))
        l = probWL(tabLettre, dicoDepSup, dicoDepInf, nbDep)
        if l > max:
            max = l
            lang = name
        print(name+"  "+str(l))
    print(filename+' appartient au '+lang+' avec un valeur bayesienne en log de '+str(max))
        
def genTexte(textLang, longN, n):
    dico = count_ngrams(textLang, n)
    dico = normaliserDict(dico)
    textG = ""
    liste = list()
    for l in dico.keys():
        dico[l] = int(dico[l]*3000)
    for l in dico.keys():
        for i in range(dico[l]):
            liste.append(l)
    random.shuffle(liste)
    #print(liste)
    for i in range(longN):
        textG += liste[random.randint(0, len(liste)-1)]
    return textG
    
def arbreB_huffman(dicoOccu):
    tas = [(cle, lettre) for lettre, cle in dicoOccu.items()]
    heapq.heapify(tas)
    while len(tas) > 1:
        cle1, lettre1 = heapq.heappop(tas)
        cle2, lettre2 = heapq.heappop(tas)
        heapq.heappush(tas, (cle1+cle2, {0:lettre1, 1:lettre2}))
    tas = tas[0][1]
    return tas

def parcoutPrefixeTasHuffman(tas, parcourt, dico_code_lettre):
    for noeud in tas:
        if len(tas[noeud]) == 1:
            dico_code_lettre[parcourt+str(noeud)] = tas[noeud]
        else:
            parcoutPrefixeTasHuffman(tas[noeud], parcourt+str(noeud), dico_code_lettre)
    return dico_code_lettre
    
def dico_huffman(tas):
    dico_code_lettre = {}
    parcoutPrefixeTasHuffman(tas,'', dico_code_lettre)
    return dico_code_lettre

#text un texte en 26 lettres d alphabet pure
def encode(text, dico_code_lettre):
    dico_lettre_code = { lettre:code for code, lettre in dico_code_lettre.items()}
    textB = ''
    for lettre in text:
        if lettre.isalpha():
            textB += dico_lettre_code[lettre]
    return textB
    
def decode(textB, dico_code_lettre):
    text = ''
    pattern = ''
    for bit in textB:
        pattern += bit
        if pattern in dico_code_lettre.keys():
            text += dico_code_lettre[pattern]
            pattern = ''
    return text
    
def compressHuffman(text):
    dicoOccu = count_ngrams(text, 1)
    arbreBH = arbreB_huffman(dicoOccu)
    #print(dicoOccu)
    #print(arbreBH)
    dico_code_lettre = dico_huffman(arbreBH)
    textB = encode(text, dico_code_lettre)
    return textB
    
def decompressHuffman(textB, dico_code_lettre):
    text = decode(textB, dico_code_lettre)
    return text
 
def compressFile(originalfile, outfilewrite):
    text = openfile(originalfile)
    text = preprocess(text)
    f = open(outfilewrite, "w")
    textB = compressHuffman(text)
    f.write(textB)
    f.close()
    
def decompressFile(originalfile, outfilewrite, dicofile):
    textB = open(originalfile, "r").read()
    f = open(outfilewrite, "w")
    fd = openfile(dicofile)
    fd = preprocess(fd)
    dicoOccu = count_ngrams(fd, 1)
    arbreBH = arbreB_huffman(dicoOccu)
    dico_code_lettre = dico_huffman(arbreBH)
    text = decompressHuffman(textB, dico_code_lettre)
    f.write(text)
    f.close
    
def esperanceCodage(dico_lettre_code, dicoProb, n):
    esperance = 0    
    for lettre in dico_lettre_code.keys():
        esperance += len(dico_lettre_code[lettre])*dicoProb[lettre]
    return n*esperance
    
def esperanceCodageText(text, n):
    dicoProb = count_ngrams(text, 1)
    dicoProb = normaliserDict(dicoProb)
    arbreBH = arbreB_huffman(dicoProb)
    dico_code_lettre = dico_huffman(arbreBH)
    dico_lettre_code = { lettre:code for code, lettre in dico_code_lettre.items()}
    return esperanceCodage(dico_lettre_code, dicoProb, n)
    
#print(divergence(p,q))
#graph()
#print(count_ngrams(open("test", "r").read(), 1))
#print(entropieText(openfile("fr3.txt"), 1))
#print(divergenceText(openfile("fr.txt"), openfile("en.txt"),1))
#histogrameFile("fr.txt",2)
#deterLangEntropie("fr4.txt",["fr3.txt", "en.txt", "es.txt", "de.txt", "it.txt"], 1)
#print(genTexte(openfile("fr.txt"), 1000 ,1))
#deterLangBaye("fr2.txt",["fr.txt", "en.txt", "es.txt","it.txt"], 4, 1)
#compressFile("fr.txt", "frB.txt")
#decompressFile("frB.txt", "frBDecomp.txt", "fr.txt")
#print(esperanceCodageText(openfile("es.txt"), 1))