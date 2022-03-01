import numpy as np 
import os
import sys
import bigfloat as big

Run = -1
Prior_Cycle = str(int(filter(str.isdigit, os.getcwd().split('/')[-1]))-1)
Prior_Run = '../Run{}'.format(Prior_Cycle)

LogLikes = np.genfromtxt('{}/True_Likelihoods.txt'.format(Prior_Run))
Likes = [big.exp(i) for i in LogLikes]
Probs = []
Sum = np.sum(Likes)
for i in Likes:
    Probs.append(float(i/Sum))
Run_Number = str(int(np.random.choice(range(1, len(Probs)+1), p = np.array(Probs))))
print(Run_Number, LogLikes[int(Run_Number)-1])
# now we jsut need to copy the output somewhere... 
import shutil
shutil.copyfile('{}/Refinement{}/refine_EM_PCS_Dual_0.pdb'.format(Prior_Run, Run_Number), '{}/Prior_Best.pdb'.format(sys.argv[1]))

Prob_Dic = {}
Prior_Data = np.genfromtxt('{}/Refinement{}/PCS_Probs_Final_Co.txt'.format(Prior_Run, Run_Number))
for i in Prior_Data:
    Prob_Dic[str(int(i[0]))] = i[2:]
outfile = open('{}/Prior_File_Co.txt'.format(sys.argv[1]), 'w')
for Res in sorted(Prob_Dic.keys(), key = lambda i:int(i)):
    outfile.write(Res+'\t'+'\t'.join([str(i) for i in Prob_Dic[Res]])+'\n')
outfile.close()

#Prob_Dic = {}
#Prior_Data = np.genfromtxt('{}/Refinement{}/PCS_Probs_Final_Fe.txt'.format(Prior_Run, Run_Number))
#for i in Prior_Data:
#    Prob_Dic[str(int(i[0]))] = i[2:]
#outfile = open('{}/Prior_File_Fe.txt'.format(sys.argv[1]), 'w')
#for Res in sorted(Prob_Dic.keys(), key = lambda i:int(i)):
#    outfile.write(Res+'\t'+'\t'.join([str(i) for i in Prob_Dic[Res]])+'\n')
#outfile.close()
