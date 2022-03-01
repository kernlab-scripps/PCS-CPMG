# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 08:30:22 2018

@author: jstiller
"""

def Crystal_Change(file_): 
    xx = open(file_, 'r') 
    lines = xx.readlines()
    xx.close()

    outfile = open(file_[:-2] + 'pdb', 'w')     
    for line in lines: 
        if 'REMARK' in line: 
            continue
        outfile.write(line) 
    outfile.close()

Crystal_Change('refine_EM_PCS_Dual_0.sa') 
