# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:58:19 2019

@author: brian
"""

import network
datapath = u"C:/users/brian/documents/github/semesterprojectce392/"
net = network.Network(datapath+"SiouxFalls_net.txt",datapath+"SiouxFalls_trips.txt")
net.ACOHeuristic()
net.userEquilibrium("FW", 200, 1e-8, net.averageExcessCost,.1)