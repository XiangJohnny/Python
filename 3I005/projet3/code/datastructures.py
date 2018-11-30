# -*- coding: utf-8 -*-
"""
Created on Fri Mar 25 16:49:17 2016

@author: Johnny
"""

import numpy as np
import pydotplus as pp
import matplotlib.pyplot as plt
import tools

class SimpleWeb:    
    def __init__(self, nbNode, name):
        """create a web graph who have nbNode node"""
        if not isinstance(nbNode, int) or nbNode < 0:
            raise tools.InvalideParamException("le parametre doit etre un entier positif")
        if not isinstance(name, str):
            raise tools.InvalideParamException("le parametre name doit etre une chaine de caractaire")
        self.nbNode = nbNode
        self.matTrans = np.zeros((nbNode, nbNode), dtype=float)
        self.tabNode = [Node(i) for i in range(nbNode)]
        self.name = name
        self.epsilonHisto = list()
        
        self.pi = None
        
        self.num = None
        self.stack = None
        self.nbFortConnexe = None
        
        self.ergodique = None
        
    def addArc(self, tail, head):
        """add a edge form head to tail in graph"""
        if (not isinstance(head, int)) or (not isinstance(tail, int)) or head < 0 or tail < 0:
            raise tools.InvalideParamException("les parametre doivent etre des entiers naturels")
        if head >= self.nbNode or tail >= self.nbNode:
            raise tools.InvalideParamException("identifiant de Noeud n est pas valide")
        arc = Arc(self.tabNode[tail], self.tabNode[head])
        if self.matTrans[tail][head] != 0:
            raise tools.InvalideParamException("arc deja present")
        self.tabNode[head].addInArc(arc)
        self.tabNode[tail].addOutArc(arc)
        
    def updateTransitionMatrix(self):
        """calcul all edge probability of transition"""
        for i in range(self.nbNode):
            node = self.tabNode[i]
            nbOutArc = len(node.tabArcOut)
            tabArc = node.tabArcOut
            poid = 0 
            if nbOutArc==0:
                self.matTrans[node.id][node.id] = 1
            else:
                poid = 1.0/nbOutArc
                for x in range(nbOutArc):
                    tabArc[x].poid=poid
                    self.matTrans[node.id][tabArc[x].head.id] = poid
                
    def getGraph(self):
        """create png format image represente this graph"""
        graph = pp.Dot(graph_type='digraph')
        for node in self.tabNode:
            graph.add_node(pp.Node(node.id))
            for arc in node.tabArcOut:
                head = arc.head.id
                tail = arc.tail.id
                graph.add_edge(pp.Edge(tail, head, label=round(self.matTrans[tail][head]-0.005, 2)))
        graph.write_png(self.name+".png")
      
    def nextStep(self, pi):
        if not isinstance(pi, list) or len(pi) != self.nbNode:
            raise tools.InvalideParamException("erreur parametre")
        pi = np.matrix(pi)
        p = np.matrix(self.matTrans)
        pp = (pi*p).tolist()
        return pp[0]
    
    def getTabEpsilon(self, pi1, pi2):
        """calcule la valeur absolue entre les deux tableaus"""
        if len(pi1) != len(pi2):
            raise tools.InvalideParamException("les tableaus sont de longeur different")
        return [abs(pi1[i]-pi2[i]) for i in range(len(pi1))]
        
    def maxEpsilon(self, tab):
        max = 0
        for i in range(len(tab)):
            if max < tab[i]:
                max = tab[i]
        return max
    
    def produitMatrixTran(self, epsilon):
        if not isinstance(epsilon, float) or epsilon <=0 or epsilon > 1:
            raise tools.InvalideParamException("erreur de parametre")
#        self.tarjan() ## dans le cas ou le graph n est pas irreductible
#        if self.nbFortConnexe != 1: #on peut dire directement que il n y a pas convergence
#            print("non convergent") #mais il y a une risque que se soit periodique
#            return #meme si il est irreductible donc regarder la variation de epsilon pour savoir s il est periodique
        self.epsilonHisto.clear()
        p = np.matrix(self.matTrans)
        p1 = np.matrix(self.matTrans)
        p0 = None
        max = 1
        limite = self.nbNode*10
        count = 0
        self.ergodique = True
        while True:
            p0 = p1
            p1 = p0*p
            ep = abs(p0-p1).max()
            if max <= ep:
                count += 1
                if count == limite:
                    self.ergodique = False
                    break
            max = ep
            self.epsilonHisto.append(max)
            if max < epsilon:
                break
        self.pi = p1.tolist()[0]
             
    def drawEpsilon(self):
        """dessiner l historique de epsilon"""
        plt.plot(self.epsilonHisto)
        plt.title(self.name+" methode pˆn\nevolotion de la convergence en fonction de nombre de pas")
        plt.xlabel("nombre de pas")
        plt.ylabel("valeur d epsilon")
        plt.show()
        
    def tarjan(self):
        
        self.nbFortConnexe = 0
        self.stack = list()
        self.num = 0
        for i in range(self.nbNode):
            self.tabNode[i].num = -1
            self.tabNode[i].numAccesible = -1
            self.tabNode[i].inStack = None
        
        def parcours(node):
            node.num = self.num
            node.numAccesible = self.num
            self.num += 1
            self.stack.append(node)
            node.inStack = True
            for i in range(len(node.tabArcOut)):
                if node.tabArcOut[i].head.num == -1:
                    parcours(node.tabArcOut[i].head)
                    node.numAccesible = min(node.numAccesible, node.tabArcOut[i].head.numAccesible)
                elif node.tabArcOut[i].head.inStack == True:
                    node.numAccesible = min(node.numAccesible, node.tabArcOut[i].head.num)
            if node.num == node.numAccesible:
                self.nbFortConnexe += 1
                while True:
                    w = self.stack.pop()
                    w.inStack = False
                    if w.id == node.id:
                        break
        for i in range(self.nbNode):
            if self.tabNode[i].num == -1:
                parcours(self.tabNode[i])
        return self.nbFortConnexe
                
class Node: 
    def __init__(self, id):
        """create node who have a integer name"""
        if not isinstance(id, int) or id < 0:
            raise tools.InvalideParamException("le parametre doit etre un entier naturel")
        self.id = id
        self.tabArcOut = list()
        self.tabArcIn = list()
        self.num = None
        self.numAccesible = None
        self.inStack = None
    
    def addOutArc(self, arc):
        """add in the list of output edge reference a new edge reference"""
        if not isinstance(arc, Arc):
            raise tools.InvalideParamException("le parametre doit etre une reference de Arc")
        self.tabArcOut.append(arc)
        
    def addInArc(self, arc):
        """add in the list of input edge reference a new edge reference"""
        if not isinstance(arc, Arc):
            raise tools.InvalideParamException("le parametre doit etre une reference de Arc")
        self.tabArcIn.append(arc)

class Arc:
    def __init__(self, tail, head):
        """create a edge"""
        if not isinstance(head, Node) or not isinstance(tail, Node):
            raise tools.InvalideParamException("head et tail doivent etre un Node")
        self.head = head
        self.tail = tail
        self.poid = None