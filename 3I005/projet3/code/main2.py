# -*- coding: utf-8 -*-
"""
Created on Sun Apr 17 01:38:00 2016
Python3
@author: Johnny
"""

import tools
import internaute as it
import simulationpi as sp
import random

if __name__ == '__main__' :
    nbNode = 500
    step = 10000
    epsilon = 0.01
    #generation d un web ergodique qvec nbNode noeuds
    nanoweb = tools.genErgodique(nbNode)
    #nanoweb.getGraph() # cree la representation image, peut prendre plusieur minute
    nanoweb.produitMatrixTran(epsilon)#Pˆn quand n->+oo 
    #nanoweb.drawEpsilon()
    #print("P^n distribution de proba pour nanoWeb1:")
    #print(nanoweb.pi)
    
    depart = random.randint(0, nbNode-1)
    bob = it.Internaute(nanoweb)
    bob.goTo(depart)
    bob.walk(step, epsilon)
    #bob.trace(1, "epsilon.txt")
    #bob.drawEpsilon()
    #print("Internaute distribution de proba pour nanoWeb1:")
    #print(bob.pi)
    
    
    sim = sp.SimulationPi(nanoweb)
    pi = [0 for i in range(nbNode)]
    pi[depart] = 1
    sim.walk(step, pi, epsilon)
    #sim.drawEpsilon()
    #print("pi*P distribution de proba pour nanoWeb1:")
    #print(sim.pi)
    #print(nanoweb.epsilonHisto)
    #print("\n")
    #print(bob.epsilonHisto)  
    #print("\n")
    #print(sim.epsilonHisto)   
    
    #comparaison entre les evolutions de epsilon
    tools.drawEpsilonCompare(nanoweb, bob, sim, 1)