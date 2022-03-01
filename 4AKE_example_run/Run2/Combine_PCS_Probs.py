import os 
import numpy as np 
import sys

Run = sys.argv[1] 
Path = 'Refinement{}/Output_Files'.format(Run)

PCS_Probs = {}

for file_ in os.listdir('{}/PCS_Calc'.format(Path)): 
    Data = np.loadtxt('{}/PCS_Calc/{}'.format(Path, file_)) 
    Res = filter(str.isdigit, file_)
    PCS_Probs[int(Res)] = [Data[-12]] 

for file_ in os.listdir('{}/Probs'.format(Path)):
    Data = np.loadtxt('{}/Probs/{}'.format(Path, file_))
    Res = filter(str.isdigit, file_)
    PCS_Probs[int(Res)].append(Data[-12])

outfile = open('Refinement{}/PCS_Probs_Final.txt'.format(Run), 'w')
outfile.write('Res\tPCS\tP1\tP2\tP3\tP4\n')
for Res in sorted(PCS_Probs.keys()):
    PC = PCS_Probs[Res][0]
    P1, P2, P3, P4 = PCS_Probs[Res][1]
    outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(Res, PC, P1, P2, P3, P4))
outfile.close() 
