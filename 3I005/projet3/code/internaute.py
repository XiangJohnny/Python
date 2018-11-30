# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 14:38:30 2016

@author: Johnny
"""
import matplotlib.pyplot as plt
import tools
import datastructures as ds
import numpy as np
import random

class Internaute:
    def __init__(self, nanoWeb):
        if not isinstance(nanoWeb, ds.SimpleWeb):
            raise tools.InvalideParamException("le parametre doit etre de type SimpleWeb")
        self.nanoWeb = nanoWeb
        self.position = None
        self.tab = [0 for i in range(nanoWeb.nbNode)]
        self.nbStep = 1
        self.epsilonHisto = list()
        self.pi = ds.Node
        
        self.ergodique = None
        
    def goTo(self, position):
        if not isinstance(position, int) or position < 0:
            raise tools.InvalideParamException("le parametre doit etre un entier naturel")
        if position >= self.nanoWeb.nbNode:
            raise tools.InvalideParamException("la destination n  existe pas")
        self.position = position
        
        
    def trace(self, step, filename):
        """enregistre dans un fichier l historique des epsilons tous les step pas"""
        if not isinstance(step, int) or step <= 0:
            raise tools.InvalideParamException("le pas doit etre un entier strictement superieur a 0")
        if not isinstance(filename, str):
            raise tools.InvalideParamException("filename doit etre une chaine de caractaire")
        if len(self.epsilonHisto) == 0:
            raise tools.InvalideAppelException("il faut appeler au moins une fois la fonction walk avant celui ci")
        f = open(filename, "w")
        s = ""
        for i in np.arange(0, len(self.epsilonHisto), step):
            s += str(i)
            s += "\t"+str(self.epsilonHisto[i])+"\n"
        f.write(s)
        f.close()
        
    def walk(self, nbStep, epsilon):
        """se deplacer dans le web au max nbStep pas et s arreter si se seuil est < epsilon"""
        if not isinstance(nbStep, int) or nbStep <= 0 or not isinstance(epsilon, float) or epsilon <=0 or epsilon > 1:
            raise tools.InvalideParamException("erreur de parametre")
        if self.position == None:
            raise tools.InvalideAppelException("il faut d abort placer l internaute  avec goTo")
        i = 0
        self.epsilonHisto.clear()
        pi1 = [0 for i in range(self.nanoWeb.nbNode)]
        pi2 = [0 for i in range(self.nanoWeb.nbNode)]
        next = -1
        max = 1
        limite = self.nanoWeb.nbNode*10
        count = 0
        self.ergodique = True
        while i < nbStep and max > epsilon:
            i += 1
            if next != -1:
                pi1[next] += 1
            next = self.choiseNext()
            self.position = next
            pi2[next] += 1
            tabEpsilon = self.nanoWeb.getTabEpsilon(self.normaliseTab(pi1), self.normaliseTab(pi2))
            ep = self.nanoWeb.maxEpsilon(tabEpsilon)
            if max <= ep:
                count += 1
                if count == limite:
                    self.ergodique = False
                    break
            max = ep
            self.epsilonHisto.append(max)
        self.nbStep = i
        self.tab = pi2
        self.pi = [ pi2[x]/i for x in range(len(pi2))]
        
    def normaliseTab(self, tab):
        sum = np.sum(tab)
        if sum == 0:
            return [0 for i in range(len(tab))]
        return [tab[i]/sum for i in range(len(tab))]
        
    def choiseNext(self):
        """etant donne notre position actuel  choisir le noeud suivant a prendre"""
        r = random.random()
        cum = 0
        tabArcOut = self.nanoWeb.tabNode[self.position].tabArcOut
        if len(tabArcOut) == 0:
            return self.position
        for i in range(len(tabArcOut)):
            cum += tabArcOut[i].poid
            if r < cum:
                return tabArcOut[i].head.id
        return tabArcOut[len(tabArcOut)-1].head.id
        
    def drawEpsilon(self):
        """dessiner l historique de epsilon"""
        plt.plot(self.epsilonHisto)
        plt.title(self.nanoWeb.name+" methode internaute\nevolotion de la convergence en fonction de nombre de pas")
        plt.xlabel("nombre de pas")
        plt.ylabel("valeur d epsilon")
        plt.show()