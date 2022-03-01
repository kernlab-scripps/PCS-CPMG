#!/bin/bash
ADK="4AKE"
cd Run1 
jid1=$(sbatch  --job-name=$ADK"_1" --parsable Run_Script.txt)
jid1b=$(sbatch -J $ADK"_1b" --dependency=afterany:$jid1 --parsable Finish_Up_Run.txt)
cd ../Run2
jid2=$(sbatch  --job-name=$ADK"_2" --parsable --dependency=afterany:$jid1b Run_Script_Follow.txt)
jid2b=$(sbatch -J $ADK"_2b" --dependency=afterany:$jid2 --parsable Finish_Up_Run.txt)
cd ../Run3
jid3=$(sbatch  --job-name=$ADK"_3" --parsable --dependency=afterany:$jid2b Run_Script_Follow.txt)
jid3b=$(sbatch -J $ADK"_3b" --dependency=afterany:$jid3 --parsable Finish_Up_Run.txt)
cd ..
sbatch -J "Finish_"$ADK --dependency=afterany:$jid3b Clean_Up_Script.txt

