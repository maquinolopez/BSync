#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 17 21:24:44 2023

@author: ma2060
"""

import PyPlum as pp                  
from pandas import read_csv           
from numpy import array, append, interp, log, hstack, mean, std, column_stack,arange, searchsorted,concatenate, load_file
from numpy.random import uniform, normal
from matplotlib.pyplot import plot, show, subplots, figure, axvline, legend
import matplotlib.lines as mlines
import pytwalk                         
from os.path import expanduser         
from os import chdir


class BSynch:
    def __init__(self, input_data, target_data, folder='~/Documents/BSynch/', 
                 thick=1., Uq = False, #model parameters 
                 mean_m=.5,shape_m=5.,mean_acc=10,shape_acc=1.5, # Autoregressive model parameters
                 iterations=2500,burnin=4000,thi=25, # twalk parameters
                 intv=.95,showchrono=False,  # plot parameters
                 Ts_mod=True, # use T distribution or normal distribution 
                 tparam=False, # tparams referes to which parammetrization to use # True: simulate alphas, False: simulates ms
                 g_thi=2,Sdate=True,seed=True,d_by=1.,plotresults=True
                 ):
        
        self.input_data = load_file(input_data, folder)
        
        
        
        
        