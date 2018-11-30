# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 14:55:36 2016
Python3
@author: Johnny
"""


import tools
import internaute as it
import simulationpi as sp
import random

if __name__ == '__main__' :
    step = 1000
    epsilon = 0.001
    
    ######################pour le nano Web 1#######################
    web1 = tools.creeNanoWeb1()
    #web1.getGraph()# cree la representation image
    web1.produitMatrixTran(epsilon)#Pˆn quand n->+oo
    web1.drawEpsilon()
    print("P^n distribution de proba pour nanoWeb1:")
    print(web1.pi)
    
    bob = it.Internaute(web1)
    #le point de depart de l internaute est aleatoire
    depart = random.randint(0, web1.nbNode-1)
    print("point de depart de nanoWeb1:"+str(depart))
    bob.goTo(depart)
    bob.walk(step, epsilon)
    #bob.trace(1, "bobepsilon.txt")
    bob.drawEpsilon()
    print("Internaute distribution de proba pour nanoWeb1:")
    print(bob.pi)
    
    sim1 = sp.SimulationPi(web1)
    #le point de depart de l internaute est aleatoire
    pi1 = [0 for i in range(web1.nbNode)]
    pi1[depart] = 1
    sim1.walk(step, pi1, epsilon)
    sim1.drawEpsilon()
    print("pi*P distribution de proba pour nanoWeb1:")
    print(sim1.pi)
    #print(web1.epsilonHisto)
    #print("\n")
    #print(bob.epsilonHisto)  
    #print("\n")
    #print(sim1.epsilonHisto)    

    #comparaison entre les evolutions de epsilon
    tools.drawEpsilonCompare(web1, bob, sim1, 1)

    ######################pour le nano Web 2#######################
    web2 = tools.creeNanoWeb2()
    web2.getGraph()# cree la representation image
    web2.produitMatrixTran(epsilon)#Pˆn quand n->+oo 
    web2.drawEpsilon()
    print("P^n distribution de proba pour nanoWeb2:")
    print(web2.pi)   

    alice = it.Internaute(web2)
    #le point de depart de l internaute est aleatoire
    depart = random.randint(0, web2.nbNode-1)
    print("point de depart de nanoWeb2:"+str(depart))
    alice.goTo(depart)
    alice.walk(step, epsilon)
    #alice.trace(1, "aliceepsilon.txt")
    alice.drawEpsilon()
    print("Internaute distribution de proba pour nanoWeb2:")
    print(alice.pi)

    sim2 = sp.SimulationPi(web2)
    #le point de depart de l internaute est aleatoire
    pi2 = [0 for i in range(web2.nbNode)]
    pi2[depart] = 1
    sim2.walk(step, pi2, epsilon)
    sim2.drawEpsilon()
    print("pi*P distribution de proba pour nanoWeb3:")
    print(sim2.pi)
    #print(web2.epsilonHisto)
    #print("\n")
    #print(alice.epsilonHisto)  
    #print("\n")
    #print(sim2.epsilonHisto)   
    
    #comparaison entre les evolutions de epsilon
    tools.drawEpsilonCompare(web2, alice, sim2, 1)    
    
    ######################pour le nano Web 3#######################
    web3 = tools.creeNanoWeb3()
    web3.getGraph()# cree la representation image
    web3.produitMatrixTran(epsilon)#Pˆn quand n->+oo 
    web3.drawEpsilon()
    print("P^n distribution de proba pour nanoWeb3:")
    print(web3.pi)

    jean = it.Internaute(web3)
    #le point de depart de l internaute est aleatoire
    depart = random.randint(0, web3.nbNode-1)
    print("point de depart de nanoWeb3:"+str(depart))
    jean.goTo(depart)
    jean.walk(step, epsilon)
    #jean.trace(1, "jeanepsilon.txt")
    jean.drawEpsilon()
    print("Internaute distribution de proba pour nanoWeb3:")
    print(jean.pi)

    sim3 = sp.SimulationPi(web3)
    #le point de depart de l internaute est aleatoire
    pi3 = [0 for i in range(web3.nbNode)]
    pi3[depart] = 1
    sim3.walk(step, pi3, epsilon)
    sim3.drawEpsilon()
    print("pi*P distribution de proba pour nanoWeb3:")
    print(sim3.pi)
    #print(web3.epsilonHisto)
    #print("\n")
    #print(jean.epsilonHisto)  
    #print("\n")
    #print(sim3.epsilonHisto)   
  
    #comparaison entre les evolutions de epsilon    
    tools.drawEpsilonCompare(web3, jean, sim3, 1)
    