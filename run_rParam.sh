#!/bin/sh

#PBS -N rParam
#PBS -l nodes=1,walltime=80:00:00
cd /home/cornejo/model/
## f value will be devided by 100, therefore if you want 
## fishing mortality of 1.2 it has to be like -f 120
##/home/cornejo/model/gridtest -lag 10 -x $i 10 -y $j 20 -f 127 -Name $k $arg
##/home/cornejo/model/gridtest -lag 0 -x 25 25 -y 50 50 -f 127 -Name rPar -Int -SST -Cur -tga
/home/cornejo/model/gridtest -lag 0 -x 25 25 -y 50 50 -f 127 -Name rPar -Int -NoSST -NoCur -tga
EOS
qstat -u cornejo
