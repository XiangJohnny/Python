# -*- coding: utf-8 -*-
"""
Created on Sat Apr 16 14:40:09 2016

@author: Johnny
"""
import tools
import datastructures as ds
import matplotlib.pyplot as plt


class SimulationPi:
    def __init__(self, nanoWeb):
        if not isinstance(nanoWeb, ds.SimpleWeb):
            raise tools.InvalideParamException("invalide parametre")
        self.nanoWeb = nanoWeb
        self.epsilonHisto = list()
        self.pi = None
        
        self.ergodique = None
    
    def walk(self, step, pi, epsilon):
        """se deplacer dans le web au max nbStep pas et s arreter si se seuil est < epsilon"""
        self.epsilonHisto.clear()
        i = 0
        pi1 = None
        pi2 = pi
        limite = self.nanoWeb.nbNode*10
        count = 0
        max = 1
        self.ergodique = True
        while i < step:
            i += 1
            pi1 = pi2
            pi2 = self.nanoWeb.nextStep(pi2)
            tabEpsilon = self.nanoWeb.getTabEpsilon(pi1, pi2)
            ep = self.nanoWeb.maxEpsilon(tabEpsilon)
            if max <= ep:
                count += 1
                if count == limite:
                    self.ergodique = False
                    break
            max = ep
            self.epsilonHisto.append(max)
            if max < epsilon:
                break
        self.pi = pi2
        
    def drawEpsilon(self):
        """dessiner l historique de epsilon"""
        plt.plot(self.epsilonHisto)
        plt.title(self.nanoWeb.name+" methode pi*P\nevolotion de la convergence en fonction de nombre de pas")
        plt.xlabel("nombre de pas")
        plt.ylabel("valeur d epsilon")
        plt.show()