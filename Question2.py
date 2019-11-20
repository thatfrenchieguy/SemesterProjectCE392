# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 17:29:50 2019

@author: brian
"""

import network
datapath = u"C:/users/brian/documents/github/semesterprojectce392/"
net = network.Network(datapath+"SiouxFalls_net.txt",datapath+"SiouxFalls_trips.txt")
Upgrades = ['(5,9)','(9,5)','(9,10)','(10,9)','(10,15)','(15,10)','(15,22)','(22,15)','(22,21)','(21,22)']
#net.ACOHeuristic()
net.userEquilibrium("FW", 200, 1e-9, net.averageExcessCost,.000001)
originalflows = dict()
for l in net.link:
    originalflows[l] = net.link[l].flow
datapath = u"C:/users/brian/documents/github/semesterprojectce392/"
net = network.Network(datapath+"SiouxFalls_net.txt",datapath+"SiouxFalls_trips.txt")
for u in Upgrades:
    net.link[u].freeFlowTime = net.link[u].freeFlowTime/2
net.userEquilibrium("FW", 200, 1e-9, net.averageExcessCost,.000001)
finalflows = dict()
for l in net.link:
    finalflows[l] = net.link[l].flow
for l in net.link:
    if finalflows[l]-originalflows[l]<-1000:
        print(l)
        print(finalflows[l]-originalflows[l])