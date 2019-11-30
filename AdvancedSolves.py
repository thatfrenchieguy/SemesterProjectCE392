# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:58:19 2019

@author: brian
"""

import network
datapath = u"C:/users/brian/documents/github/semesterprojectce392/"
net = network.Network(datapath+"SiouxFalls_net.txt",datapath+"SiouxFalls_trips.txt")
#hot start
net.ACOHeuristic()
net.ACOWriteOut()
#FW, max iterations, AEC gap tolerance, AEC as choice of method, Frank Wolve precision
net.userEquilibrium("FW", 20000, 1e-5, net.averageExcessCost,.0000000001)