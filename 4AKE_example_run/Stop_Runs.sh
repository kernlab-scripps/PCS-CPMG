ADK="4AKE"
for i in `seq 1 5`;
do
    scancel --name="JBS_"$ADK"_"$i
    scancel --name="JBS_"$ADK"_"$i"b"
done