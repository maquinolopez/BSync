#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 13:52:45 2023

@author: ma2060
"""

import numpy as np 
from scipy.interpolate import interp1d
from scipy.stats import t, gamma, beta
from pytwalk import pytwalk
import matplotlib.pyplot as plt

def BSynch(input_data, target_data, folder='~/Documents/BSynch/', 
           ta=3, tb=4,  
           shape_acc=10, meanM=0.7, alpha=4, 
           iters=2000, burn=1000, thin=150):

    # Load input data
    inp = load_file(input_data, folder)
    
    depth = inp['MCD']
    age_temp = inp['Age'] * 1000
    
    if 'ProxyType' in inp:
        # Handle proxy transformations
        if inp['ProxyType'][0] in ['d18op', 'd18ob']:
            inp['ProxyValue'] = -inp['ProxyValue'] 
        elif inp['ProxyType'][0] == 'mg':
            inp['ProxyValue'] = np.log(inp['ProxyValue'])
            
    inp = inp[inp['Age'] <= max(target_data['Age']) - 10]
    
    # Load target data
    tar = load_file(target_data, folder)
    
    # Scale input and target data
    inp['ProxyValue'] = (inp['ProxyValue'] - inp['ProxyValue'].min()) / inp['ProxyValue'].ptp()
    tar['ProxyValue'] = (tar['ProxyValue'] - tar['ProxyValue'].min()) / tar['ProxyValue'].ptp()
    
    # Set up model parameters
    N = int(np.round((max(inp['Age']) - min(inp['Age']))/100)) 
    interval = (max(inp['Age']) - min(inp['Age'])) / N
    nodes = np.linspace(min(inp['Age']), max(inp['Age']), N+1)
    sr = abs((nodes[-1] - nodes[0]) / (max(tar['Age']) - min(tar['Age'])))
    scale_acc = meanM / shape_acc

    # Define priors 
    shape_acc = shape_acc
    scale_acc = meanM / shape_acc
    beta = (alpha / meanM) - alpha
    
    # Define likelihood
    def likelihood(params):
        # Accumulation model
        d = np.zeros(N)
        d[0] = params[2]
        for i in range(1, N):
            d[i] = params[1] * d[i-1] + (1 - params[1]) * params[i+2]
        
        points = np.concatenate([[params[0]], params[0] + np.cumsum(interval / d)]) 
        f = interp1d(nodes, points)
        new_input = np.c_[f(inp['Age']), inp['ProxyValue']]
        ft = interp1d(tar['Age'], tar['ProxyValue'])
        new_target = ft(new_input[:,0])
        
        # Likelihood
        resid = new_input[:,1] - new_target
        return np.sum(t.logpdf(resid, df=tb, loc=0, scale=sigma))

    # Prior density
    def prior_density(params):
        dstart = t.logpdf(params[0], df=2, loc=min(tar['Age']), scale=5)
        dmem = beta.logpdf(params[1], a=alpha, b=beta)
        drate = gamma.logpdf(params[2:], a=shape_acc, scale=scale_acc)
        return dstart + dmem + np.sum(drate)

    # Objective
    def objective(params):
        return -(prior_density(params) + likelihood(params))
    
    # Proposal distribution
    def propose():
        # Sample params
        return params
    
    # Run sampler
    output = pytwalk(objective, propose, iters, burn, thin) 
    samples = output[burn::thin]

    # Extract results
    accrates = np.zeros((len(samples), N-1))
    for i,sample in enumerate(samples):
       accrates[i,:] = sample[3:] 

    ages = np.zeros((len(samples), N+1)) 
    for i,sample in enumerate(samples):
        ages[i,:] = np.concatenate([[sample[0]],  
                                    sample[0] + np.cumsum(interval / accrates[i,:])])

    memory = np.quantile(samples[:,1], [0.05, 0.5, 0.95])
    quants = np.quantile(ages, [0.05, 0.32, 0.5, 0.68, 0.95], axis=0)

    # Synchronized input
    align = np.c_[interp1d(nodes, quants[1,:])(inp['Age']), inp['ProxyValue']]

    # Plots
    fig, axs = plt.subplots(2, 2, figsize=(10, 8))

    # Prior vs posterior plots
    axs[0,0].plot(gamma.pdf(x, a=shape_acc, scale=scale_acc)) 
    axs[0,0].plot(bkde(accrates[:,:-1].flatten()))
    axs[0,1].plot(beta.pdf(x, a=alpha, b=beta))
    axs[0,1].plot(bkde(samples[:,1]))

    # Age-depth  
    axs[1,0].plot(quants[1,:], interp1d(age_temp, depth)(quants[1,:]))
    axs[1,0].fill_between(quants[0,:], quants[4,:], alpha=0.5)

    # Aligned data
    axs[1,1].plot(align[:,0], align[:,1])

    fig.tight_layout()
    plt.show()

    return samples, ages, memory

# Load data
def load_file(filename, folder):
    path = os.path.join(folder, filename + '.csv')
    if os.path.exists(path):
        return pd.read_csv(path)
    path = os.path.join(folder, filename + '.txt') 
    if os.path.exists(path):
        return pd.read_csv(path, sep='\t')
    print('No file found for {}!'.format(filename))
    return None

# Kernel density estimation
from sklearn.neighbors import KernelDensity
def bkde(x):
    kde = KernelDensity(bandwidth=0.5).fit(x.reshape(-1, 1))
    xgrid = np.linspace(x.min(), x.max(), 1000)[:,None]
    logprob = kde.score_samples(xgrid)
    return xgrid, np.exp(logprob)

# Sample from proposal
def propose():
    dstart = np.random.uniform(min(tar['Age']), min(tar['Age']) + 10) 
    mem = np.random.uniform(0.5, 1)
    arate = np.random.uniform(0.5*sr, 1.5*sr, size=N-1)
    return np.concatenate([[dstart], [mem], arate])

# Plotting functions
import matplotlib.pyplot as plt
def traceplot(samples):
    fig, axs = plt.subplots(ncols=samples.shape[1], sharey=True, figsize=(12, 4))
    for i in range(samples.shape[1]):
        axs[i].plot(samples[:,i])
        axs[i].set_title(f'Parameter {i+1}')
    fig.tight_layout()

def pairplot(samples):
   pd.DataFrame(samples).plot.scatter_matrix(diagonal='hist')