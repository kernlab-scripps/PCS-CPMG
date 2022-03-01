import os 
Ref_List = []
for Entry in os.listdir('.'):
    if 'Refinement' in Entry:
        Ref_List.append(int(filter(str.isdigit, Entry)))
Run_Array = sorted(Ref_List)

Stats, Likes = [], []
for Run in Run_Array:
    try:
        with open('Refinement{}/Statistics.txt'.format(Run), 'r') as fh: 
            if Stats == []:
                Stats.append(''.join(fh.readlines()))
            else:
                Stats.append(fh.readlines()[2])
        with open('Refinement{}/True_Likelihoods.txt'.format(Run), 'r') as fh: 
            Likes.append(fh.readlines()[0])
    except:
	with open('Failures.txt', 'a+') as out: 
	    out.write('Problem with Run {}'.format(Run))

with open('Statistics.txt', 'w') as out: 
    out.write(''.join(Stats))

with open('True_Likelihoods.txt', 'w') as out: 
    out.write('\n'.join(Likes))
