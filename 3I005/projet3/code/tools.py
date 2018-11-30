# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 17:23:30 2016

@author: Johnny
"""

import random
import numpy as np
import datastructures as ds
import matplotlib.pyplot as plt

def creeNanoWeb1():
    n = ds.SimpleWeb(10, "nano_web1") # 10 noeuds de 0 a 9
    n.addArc(0, 1); n.addArc(0 ,4)
    n.addArc(1, 2)
    n.addArc(2, 3); n.addArc(2, 4)
    n.addArc(3, 9)
    n.addArc(4, 2); n.addArc(4, 5); n.addArc(4, 6)
    n.addArc(5, 6)
    n.addArc(6, 5); n.addArc(6, 7)
    n.addArc(7, 8)
    n.addArc(8, 7)  
    n.updateTransitionMatrix()
    return n
    
def creeNanoWeb2():
    n = ds.SimpleWeb(10, "nano_web2")
    n.addArc(0, 9)
    n.addArc(1, 0); n.addArc(1, 5)
    n.addArc(2, 1); n.addArc(2, 4)
    n.addArc(3, 2)
    n.addArc(4, 3)
    n.addArc(5, 4)
    n.addArc(6, 5)
    n.addArc(7, 6)
    n.addArc(8, 7)
    n.addArc(9, 8); n.addArc(9, 2)
    n.updateTransitionMatrix()
    return n
    
def creeNanoWeb3():
    n = ds.SimpleWeb(10, "nano_web3")
    n.addArc(0, 1)
    n.addArc(1, 2); n.addArc(1, 3)
    n.addArc(2, 3); n.addArc(2, 9)
    n.addArc(4, 5)
    n.addArc(5, 4)
    n.addArc(6, 7)
    n.addArc(7, 8)
    n.addArc(8, 7)
    n.updateTransitionMatrix()
    return n
    
def drawEpsilonCompare(nanoweb,internaute, sim, step):
    eh1 = nanoweb.epsilonHisto
    eh2 = internaute.epsilonHisto
    eh3 = sim.epsilonHisto
    x = [x for x in np.arange(0, max(len(eh1), len(eh2), len(eh3)), step)]
    
    plt.ylim(0, 1.1)
    y1 = [eh1[i] for i in np.arange(0, len(eh1), step)]
    y2 = [eh2[i] for i in np.arange(0, len(eh2), step)]
    y3 = [eh3[i] for i in np.arange(0, len(eh3), step)]
    if nanoweb.ergodique:
        labnanoWeb = "puissance P, Pas: "+ str(len(eh1))        
    else:
        labnanoWeb = "puissance P, non convergent"
    if internaute.ergodique:
        labInt = "Internaute, Pas: "+ str(len(eh2))
    else:
        labInt = "Internaure, non congergent"
    if sim.ergodique:
        labSim = "Pi*P, Pas: "+str(len(eh3))
    else:
        labSim = "Pi*P, non convergent"
    plt.plot(x[:len(y1)], y1, color="blue", label=labnanoWeb)
    plt.plot(x[:len(y2)], y2, color="red", label=labInt)
    plt.plot(x[:len(y3)], y3, color="black", label=labSim)
    plt.legend(loc='upper right', frameon=False)
    plt.title(str(nanoweb.nbNode)+" Nodes\n Evolotion de la convergence\nen fonction de nombre de pas\navec echantionnage pour tous les "+str(step)+" pas")
    plt.xlabel("nombre de pas")
    plt.ylabel("valeur d epsilon")
    plt.show()
    
def genErgodique(nbNode):
    """generateur d une nanoWeb ergodique de nbNode noeuds"""
    if not isinstance(nbNode, int) or nbNode < 3:
        raise InvalideParamException("nbNode doit etre >= 3")
    n = ds.SimpleWeb(nbNode, "ergodique graph")
    ensemble = [[n.tabNode[i]] for i in range(n.nbNode)]
    #apediodique avec construction de periode 2 et 3
    e1 = ensemble[0]
    n1 = e1[0]
    e2 = ensemble.pop(1)
    e3 = ensemble.pop(1)
    n.addArc(n1.id, e2[0].id)
    n.addArc(e2[0].id, e3[0].id)
    n.addArc(e3[0].id, n1.id)
    
    n.addArc(e2[0].id, n1.id)
    e1.extend(e2)
    e1.extend(e3)
    
    i1 = None
    i2 = None
    while len(ensemble) > 1:
        i1 = random.randint(0, len(ensemble)-1)
        while True:
            i2 = random.randint(0, len(ensemble)-1)
            if i1 != i2:
                break
        e1 = ensemble[i1]
        e2 = ensemble.pop(i2)
        n1 = e1[random.randint(0, len(e1)-1)]
        n2 = e2[random.randint(0, len(e2)-1)]
        n.addArc(n1.id, n2.id)
        n1 = e1[random.randint(0, len(e1)-1)]
        n2 = e2[random.randint(0, len(e2)-1)]
        n.addArc(n2.id, n1.id)
        e1.extend(e2)
    n.updateTransitionMatrix()
    return n
        
        
class InvalideParamException(Exception):
    def __init__(self, msg):
        self.message = msg
        
class InvalideAppelException(Exception):
    def __init__(self, msg):
        self.message = msg