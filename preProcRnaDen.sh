#!/bin/bash
rm data.csv data.temp data.temp.1 data.temp.2 DES.out Energy.out Energy.tempy RMS.out RG.out RG.temp
rnaDenOutFile='rna.denovo.1HC8.native.helix.test.650.spectime.out'
outF='rnaden.test.native.helix.spectime.650.out'
grep "SCORE" $rnaDenOutFile > score.out
scoreF='score.out'
awk -F '  ' '{print $2}' $scoreF > Energy.tempy
sed '1d' Energy.tempy > Energy.out
sed -i '1iEnergy' Energy.out
awk -F ' ' '{print $16}' $scoreF > RMS.out
awk -F ' ' '{print $28}' $scoreF > DES.out
grep "rna_rg" $outF > RG.temp
sed 's/^..............................//g' RG.temp | awk -F '      ' '{print $2}' > RG.out
sed -i "1iRG" RG.out
:|paste -d', ' DES.out - Energy.out > data.temp
:|paste -d', ' data.temp - RMS.out > data.temp.1
:|paste -d', ' data.temp.1 - RG.out > data.temp.2
awk -F '' '{if($2 > 20) {print}}' data.temp.2 | head -1
sed 's/\s*,\s*/,/g' data.temp.2 > data.csv
head data.csv