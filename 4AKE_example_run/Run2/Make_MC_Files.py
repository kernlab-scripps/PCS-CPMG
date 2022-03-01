# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 14:29:08 2018

@author: jstiller
"""

#==============================================================================
# Modules and Flags
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os

font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 15}
matplotlib.rc('font', **font)
axis_font = {'fontname':'Arial', 'weight' : 'bold', 'size':'15'}
title_font = {'fontname':'Arial', 'weight' : 'bold', 'size':'18'}
text_font = {'fontname':'Arial', 'weight' : 'bold', 'size':'15'}

#==============================================================================
# Modules and Flags
#==============================================================================
def Adjust_Prob(Square_Dic, Sigma):
    New_Prob_Dic = {}
    for Res in sorted(Square_Dic.keys(), key= lambda i:int(i)):
        Prob_List = []
        s1 = 1.0/(Sigma*np.sqrt(2*np.pi))
        s2 = -1.0/(2.0*np.power(Sigma,2))
        for Ent, Val in enumerate(Square_Dic[Res]):
            Squared_Diff = Square_Dic[Res][Ent]
            Prob_List.append(s1*np.exp(s2*Squared_Diff))
        New_Prob_Dic[Res] = [float(i)/sum(Prob_List) for i in Prob_List]
    return New_Prob_Dic
def Calc_Sigma(Square_Dic, Prob_Dic):
    Sigma_All = np.array([])
    for Res in sorted(Square_Dic.keys(), key = lambda i:int (i)):
        Indiv_Sum = np.array([])
        for Ent, Val in enumerate(Square_Dic[Res]):
            Squared_Diff = Square_Dic[Res][Ent]*Prob_Dic[Res][Ent]
            Indiv_Sum = np.append(Indiv_Sum, Squared_Diff)
        Sigma_All = np.append(Sigma_All, np.sum(Indiv_Sum))
    return np.sqrt(np.sum(Sigma_All)/(len(Sigma_All)-1))
def Calculate_Class_Av(Prob_Dic):
    Averages_List = []
    Temp_Dic = {}
    for Res in Prob_Dic:
        for ent, val in enumerate(Prob_Dic[Res]):
            Temp_Dic.setdefault(ent, []).append(val)
    Averages_List = [np.average(np.array(Temp_Dic[i])) for i in Temp_Dic]
    for Res in Prob_Dic:
        Temp_Probs = np.array(Prob_Dic[Res])*np.array(Averages_List)
        Prob_Dic[Res] = [float(i)/sum(Temp_Probs) for i in Temp_Probs]
    return Prob_Dic
def Determine_Square_Dic(Run):
    Square_Dic = {}
    fname = 'Refinement{}/refine_EM_PCS_0.sa.viols'.format(Run)
    with open(fname, 'rb') as fh:
        lines = fh.readlines()
        for ent,line in enumerate(lines):
            if 'ZN' in line and 'True' in line:
                s = filter(None, line.split(' ')) 
                Square_Dic.setdefault(s[-8],[]).append(float(s[-3])**2)
    return Square_Dic
def Calc_Total_Like(Square_Dic, Sigma, Prob_Dic):
    Q_All = np.array([])
    s1 = (-0.5*np.log(2*np.pi))
    s2 = (-1.0)*np.log(Sigma)
    s3 = (-0.5)*(1/np.power(Sigma,2))
    for Res in sorted(Square_Dic.keys(), key = lambda i: int(i)):
        Indiv_Sum = np.array([])
        for Ent, Val in enumerate(Square_Dic[Res]):
            Squared_Diff = Square_Dic[Res][Ent]
            Indiv_Sum = np.append(Indiv_Sum, Prob_Dic[Res][Ent]*(s1+s2+s3*Squared_Diff))
        Q_All = np.append(Q_All, np.sum(Indiv_Sum))
    return np.sum(Q_All)
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
    Class_Average = np.mean([np.array(Prob_Dic[i]) for i in Prob_Dic], axis = 0)
    Q_All = []
    s1 = (-0.5*np.log(2*np.pi))
    s2 = (-1.0)*np.log(Sigma)
    s3 = (-0.5)*(1/np.power(Sigma,2))
    for Res in sorted(Square_Dic.keys(), key = lambda i: int(i)):
        Indiv_Sum = np.array([])
        for Ent, Val in enumerate(Square_Dic[Res]):
            Squared_Diff = Square_Dic[Res][Ent]
            Indiv_Sum = np.append(Indiv_Sum, (s1+s2+s3*Squared_Diff+np.log(Class_Average[Ent])))
        Q_All.append(big.log(np.sum([big.exp(i) for i in Indiv_Sum])))
    return np.sum(Q_All)
    
    
def Determine_Final_Sigma(Run):
    fname = 'Refinement{}/QS_Output.txt'.format(Run)
    with open(fname, 'rb') as fh:
        return float(fh.readlines()[-12].split('\t')[1])

def Calc_Final_Q(Run):
    Q_List, Sigma_S = [-1000], Determine_Final_Sigma(Run)
    Sigma_L = []
    Square_Dic = Determine_Square_Dic(Run)
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
        

Ref_List = []
for Entry in os.listdir('.'):
    if 'Refinement' in Entry: 
        Ref_List.append(int(filter(str.isdigit, Entry)))
        
Q_Final = []
for Run in sorted(Ref_List):
    Q_True = Calc_Final_Q(str(Run))
    Q_Final.append(Q_True)

Likelihoods = [big.exp(i) for i in Q_Final]
Norm_Like = [i/np.sum(Likelihoods) for i in Likelihoods]

outfile = open('True_Likelihoods.txt', 'w')
for ent, Run in enumerate(sorted(Ref_List)):
    outfile.write('{}\t{}\t{}\n'.format(Run, Q_Final[ent], Norm_Like[ent]))
outfile.close()
