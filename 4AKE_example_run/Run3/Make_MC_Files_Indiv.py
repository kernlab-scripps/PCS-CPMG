# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:29:08 2018

@author: jstiller
"""

#==============================================================================
# Modules and Flags
#==============================================================================
import numpy as np
import math

Metals = ['Co', 'Fe']

#==============================================================================
# Modules and Flags
#==============================================================================
def Adjust_Prob(Square_Dic, Sigma):
    New_Prob_Dic = {M: {} for M in Metals}
    s1 = 1.0/(Sigma*np.sqrt(2*np.pi))
    s2 = -1.0/(2.0*np.power(Sigma, 2))
    for M in Metals: 
        for Res in Square_Dic[M]:
            Prob_List = s1*np.exp(s2*np.array(Square_Dic[M][Res])) + 1e-50
            New_Prob_Dic[M][Res] = [float(i)/sum(Prob_List) for i in Prob_List]
    return New_Prob_Dic


def Calc_Sigma(Square_Dic, Prob_Dic):
    Sigma = []
    for M in Metals: 
        for Res in Square_Dic[M]:
            Sigma.append(np.sum(np.array(Square_Dic[M][Res])*Prob_Dic[M][Res]))
    Sigma_Sum = np.sum(Sigma)
    return math.sqrt(Sigma_Sum / (len(Sigma)-1.0))


def Calculate_Class_Av(Prob_Dic):
    # Determine Class Averages:
    Prob_Array = []
    for M in Metals:
        for Res in Prob_Dic[M]:
            Prob_Array.append(Prob_Dic[M][Res])
    Class_Averages = np.mean(np.array(Prob_Array), axis = 0)
    # Adjust Probability Dic with class averages and normalize values:
    for M in Prob_Dic:
        for Res in Prob_Dic[M]:
            Prob_Dic[M][Res] *= Class_Averages
            Prob_Dic[M][Res] /= np.sum(Prob_Dic[M][Res])
    return Prob_Dic


def Determine_Square_Dic():
    Square_Dic = {M: {} for M in Metals}
    fname = 'refine_EM_PCS_Dual_0.sa.viols'
    with open(fname, 'r') as fh:
        lines = fh.readlines()
        for ent, line in enumerate(lines):
            if 'ZN' in line and 'True' in line:
                s = list(filter(None, line.split(' ')))
                s2 = list(filter(None, lines[ent-4].split(' ')))[-4]
                M = s2.split('_')[0]
                Square_Dic[M].setdefault(s[-8],[]).append(float(s[-3])**2)
    return Square_Dic


def Calc_Total_Like(Square_Dic, Sigma, Prob_Dic):
    # Calculate Total Probability Weighted Log Likelihood - Q
    Q_All = 0
    # Prepare Log Gaussian Likelihood Terms
    s1 = -0.5*np.log(2.0*np.pi)
    s2 = -1.0*np.log(Sigma)
    s3 = -0.5/np.square(Sigma)
    # Calculate Log Likelihood on Individaul Residue and add to Q
    for M in Metals: 
        for Res in Square_Dic[M]:
            Squared_Array = np.array(Square_Dic[M][Res])
            Q_All += np.sum((s1+s2+s3*Squared_Array) * Prob_Dic[M][Res])
    return Q_All


def Calc_Total_Like_Correct(Square_Dic, Sigma, Prob_Dic):
    Class_Average = np.mean([np.array(Prob_Dic[i]) for i in Prob_Dic], axis = 0)
    Q_All = np.array([])
    s1 = (-0.5*np.log(2*np.pi))
    s2 = (-1.0)*np.log(Sigma)
    s3 = (-0.5)*(1/np.power(Sigma,2))
    for Res in sorted(Square_Dic.keys(), key = lambda i: int(i)):
        Indiv_Sum = np.array([])
        for Ent, Val in enumerate(Square_Dic[Res]):
            Squared_Diff = Square_Dic[Res][Ent]
            Indiv_Sum = np.append(Indiv_Sum, Prob_Dic[Res][Ent]*(s1+s2+s3*Squared_Diff+np.log(Class_Average[Ent])))
        Q_All = np.append(Q_All, np.sum(Indiv_Sum))
    return np.sum(Q_All)
    

import bigfloat as big
def Calc_True_Likelihood(Square_Dic, Sigma, Prob_Dic):
    Prob_Array = []
    for M in Metals:
        for Res in Prob_Dic[M]:
            Prob_Array.append(Prob_Dic[M][Res])
    Class_Averages = np.mean(np.array(Prob_Array), axis = 0)
    Q_All = []
    s1 = (-0.5*np.log(2*np.pi))
    s2 = (-1.0)*np.log(Sigma)
    s3 = (-0.5)*(1/np.power(Sigma,2))
    for M in Metals: 
        for Res in Square_Dic[M]:
            Indiv_Sum = (s1+s2+s3*np.array(Square_Dic[M][Res]) + 
                         np.log(Class_Averages))
            Q_All.append(big.log(np.sum([big.exp(i) for i in Indiv_Sum])))
    return np.sum(Q_All)
    
    
def Determine_Final_Sigma():
    fname = 'QS_Output.txt'
    with open(fname, 'r') as fh:
        return float(fh.readlines()[-1].split('\t')[1])


def Calc_Final_Q():
    Q_List, Sigma_S = [-1000], Determine_Final_Sigma()
    Sigma_L = []
    Square_Dic = Determine_Square_Dic()
    Prob_Dic = Adjust_Prob(Square_Dic, Sigma_S)
    Sigma = Calc_Sigma(Square_Dic, Prob_Dic)
    Q_List.append(Calc_Total_Like(Square_Dic, Sigma, Prob_Dic))
    while abs(Q_List[-1]/Q_List[-2] - 1) > 0.0001:
        Prob_Dic = Calculate_Class_Av(Prob_Dic)
        Sigma = Calc_Sigma(Square_Dic, Prob_Dic)
        Prob_Dic = Adjust_Prob(Square_Dic, Sigma)
        Prob_Dic = Calculate_Class_Av(Prob_Dic)
        Sigma = Calc_Sigma(Square_Dic, Prob_Dic)
        Prob_Dic = Adjust_Prob(Square_Dic, Sigma)
        Q_List.append(Calc_True_Likelihood(Square_Dic, Sigma, Prob_Dic))
        Sigma_L.append(Sigma)
    return Q_List[-1]
        
Q_True_Log = Calc_Final_Q()
#Q_True = big.exp(Q_True_Log)

outfile = open('True_Likelihoods.txt', 'w')
outfile.write('{}'.format(Q_True_Log))
outfile.close()
