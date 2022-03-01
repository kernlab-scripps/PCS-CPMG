mkdir -p "Stat_LikeDir"
for i in `seq 1 5`
do 
    cp "Run"$i"/Statistics.txt" "Stat_LikeDir/Statistics"$i".txt"
    cp "Run"$i"/True_Likelihoods.txt" "Stat_LikeDir/True_Likelihoods"$i".txt"
done

