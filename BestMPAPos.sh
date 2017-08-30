#!/bin/bash

for ((a = 15; a <= 30; a = a+5))
do
  for ((b = 15; b <= 90; b = b+5))
  do
    vary="$a-$b"
    cat <<EOS | qsub -
#!/bin/sh
#PBS -N mpa_$vary
#PBS -l nodes=1
#PBS -l walltime=70:00:00

cd /home/cornejo/model/
/home/cornejo/model/gridtest -x $a 5 -y $b 5 -Int -SST -Cur -F 122
EOS
  done
done               


exit 0
