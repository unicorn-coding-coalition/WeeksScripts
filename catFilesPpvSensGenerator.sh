declare -a fTypes=minCorr10.minCount5 
minCorr5.minCount5 
win3.minCorr10.minCount5
win3.minCorr5.minCount5
win4.minCorr10.minCount5
win4.minCorr5.minCount5 )
for idx in ${!fTypes[*]};
	do
	for i in $(ls | grep "${fTypes[$idx]}");
		do 
			cat $i >> ${fTypes[$idx]}.ps.mean;
		done
	done
	