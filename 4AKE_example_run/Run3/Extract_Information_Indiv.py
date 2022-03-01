"""
Created on Wed May 16 18:18:14 2018

@author: jstiller
"""

#==============================================================================
# Goals
#==============================================================================
# Goals: First Extract RMSD
# Second: Organize Entries via RMSD
# Third: Caluclate how on Good the PCS fit is 
# Relate how good the PCS fit is with how 'correct it is' 

import os 
ADK="4AKE"

#==============================================================================
# Modules
#==============================================================================
import numpy as np

#==============================================================================
# Extract and Organize RMSD
#==============================================================================
def Extract_RMSD(): 
    file_ = open('Thesius_Alignment/theseus_sup.pdb', 'r')
    lines = file_.readlines()
    file_.close()
    for line in lines:
        if 'Classical pairwise' in line and 'REMARK' in line: 
            return float(filter(None, line.split(' '))[-1])
   
RMSD = Extract_RMSD() 

#==============================================================================
# Metal Difference
#==============================================================================
def Extract_Metal_Coor(file_):
    xx = open(file_, 'r')
    lines = xx.readlines()
    xx.close()
    
    N_Coor = {}
    for line in lines: 
        if 'ZN' in line and 'ATOM' in line:
            if 'refine' in file_:
                Metal_Coor = [float(i) for i in filter(None, line.split(' '))[5:8]]
            else:
                Metal_Coor = [float(i) for i in filter(None, line.split(' '))[6:9]]
        split_line = filter(None, line.split(' '))
        try:
            if split_line[2] == 'N':
                if 'XPLOR' in file_:
                    N_Coor.setdefault(split_line[4], np.array([float(i) for i in filter(None, line.split(' '))[5:8]]))
                else:
                    N_Coor.setdefault(split_line[5], np.array([float(i) for i in filter(None, line.split(' '))[6:9]]))
        except:
            pass
    return Metal_Coor, N_Coor
        
def Calculate_Distance(Metal1, Metal2):
    return np.sqrt(np.sum(np.power(np.array(Metal1) - np.array(Metal2),2)))


M1, N1 = Extract_Metal_Coor('Thesius_Alignment/theseus_refine_EM_PCS_Dual_0.pdb')
M2, N2 = Extract_Metal_Coor('Thesius_Alignment/theseus_{}_Modeled.pdb'.format(ADK))
RMS_M = Calculate_Distance(M1, M2)


#==============================================================================
# ATP and AMP Lid Positions
#==============================================================================  
def Calc_Angle(P1, P2, Ori):
    v1 = P1 - Ori
    v2 = P2 - Ori
    Mag1 = np.sqrt(np.sum(np.power(P1 - Ori,2)))
    Mag2 = np.sqrt(np.sum(np.power(P2 - Ori,2)))
    DotPro = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    return np.arccos(DotPro / (Mag1*Mag2))*180/np.pi

def Calc_Dist(Structure):
    # For Center of Mass Distance Measurements
    AMP_Lid = [str(i) for i in range(30, 54)]
    ATP_Lid = [str(i) for i in range(127, 160)]
    Core_Domain = range(2,31) + range(63, 97) +range(107,123) + range(160,211)
    Core_Domain = [str(i) for i in Core_Domain]
    Core, AMP, ATP = [],[],[]

    # For Angle Determintion
    Core_A = range(1,9)  + range(79,86) + range(108,115) + range(190,199)
    Core_A = [str(i) for i in Core_A]
    Lid_A = [str(i) for i in range(127,160)]
    Hinge_A = [str(i) for i in range(165,170)]
    NMP_A = [str(i) for i in range(50,60)]
    Core_R, Lid_R, Hinge_R, NMP_R = [],[],[],[]


    if 'refine' in Structure: 
        Res_Pos = 4
        Array = [5,8]
        line_count = 12
    else: 
        Res_Pos = 5
        Array = [6,9]
        line_count = 13

    with open(Structure, 'r') as fh:
        for line in fh: 
            s = filter(None, line.split(' '))
            if len(s) != line_count:
                continue
            if s[2] not in ['N', 'CA', 'C', 'O']:
                continue
            if s[Res_Pos] in AMP_Lid and s[0] == 'ATOM': 
                AMP.append(np.array([float(i) for i in s[Array[0]:Array[1]]]))
            elif s[Res_Pos] in ATP_Lid and s[0] == 'ATOM': 
                ATP.append(np.array([float(i) for i in s[Array[0]:Array[1]]]))
            elif s[Res_Pos] in Core_Domain and s[0] == 'ATOM':
                Core.append(np.array([float(i) for i in s[Array[0]:Array[1]]]))
            
            if s[Res_Pos] in Core_A and s[0] == 'ATOM': 
                Core_R.append(np.array([float(i) for i in s[Array[0]:Array[1]]])) 
            elif s[Res_Pos] in Lid_A and s[0] == 'ATOM': 
                Lid_R.append(np.array([float(i) for i in s[Array[0]:Array[1]]])) 
            elif s[Res_Pos] in Hinge_A and s[0] == 'ATOM': 
                Hinge_R.append(np.array([float(i) for i in s[Array[0]:Array[1]]])) 
            elif s[Res_Pos] in NMP_A and s[0] == 'ATOM': 
                NMP_R.append(np.array([float(i) for i in s[Array[0]:Array[1]]])) 
                
                

    Core_C = np.array(np.matrix(np.array(Core)).mean(0))[0]
    AMP_C = np.array(np.matrix(np.array(AMP)).mean(0))[0]
    ATP_C = np.array(np.matrix(np.array(ATP)).mean(0))[0]

    C_M = np.sqrt(np.sum(np.power(AMP_C - Core_C,2)))
    C_T = np.sqrt(np.sum(np.power(ATP_C - Core_C,2)))

    Core_RC = np.array(np.matrix(np.array(Core_R)).mean(0))[0]
    Lid_RC = np.array(np.matrix(np.array(Lid_R)).mean(0))[0]
    Hinge_RC = np.array(np.matrix(np.array(Hinge_R)).mean(0))[0]
    NMR_RC = np.array(np.matrix(np.array(NMP_R)).mean(0))[0]

    Theta1 = Calc_Angle(Core_RC,Lid_RC, Hinge_RC)
    Theta2 = Calc_Angle(NMR_RC, Hinge_RC, Core_RC)

    return C_M, C_T, Theta1, Theta2

file_ = 'refine_EM_PCS_Dual_0.pdb'
AMP, ATP, Theta1, Theta2 = Calc_Dist(file_)

y,x, T1, T2= Calc_Dist('{}_Modeled.pdb'.format(ADK))

outfile = open('Statistics.txt', 'w')
outfile.write('RMSD\tRMS_M\tAMP\tATP\tTheta1\tTheta2\n')
outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('NaN', 'NaN', y, x, T1, T2))
outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(RMSD, 
                  RMS_M, AMP, ATP, Theta1, Theta2))
outfile.close()

    
    
