#!/bin/sh

if [ -z $1 ]
then
	echo " Int, SST and Cur not Defined!!! "
	echo " Nothing to do "
	exit 0
fi

echo "$1 $2 $3"

if [ "$1" == "Yes" -o "$1" == "yes" ]
then
	if [ "$2" == "Yes" -o "$2" == "yes" ]
	then
		if [ "$3" == "Yes" -o "$3" == "yes" ]
		then
			walltime="60:00:00"
			arg="-Int -SST -Cur"
			echo "1 $arg"
		else
			walltime="50:00:00"
			arg="-Int -SST -NoCur"
			echo "2 $arg"
		fi
	else
		walltime="40:00:00"
		arg="-Int -NoSST -NoCur"	
		echo "3 $arg"
	fi
else
	walltime="20:00:00"
	arg="-NoSST -NoCur -NoInt"
	echo "4 $arg"
fi

cX=("20")
cY=("40")  	
for i in ${cX[@]}
do
	for j in ${cY[@]}
	do
		for k in {34..34}
		do

cat <<EOS | qsub -
#!/bin/sh
#PBS -N $k
#PBS -l nodes=1,walltime=$walltime
cd /home/cornejo/model/
## f value will be devided by 100, therefore if you want 
## fishing mortality of 1.2 it has to be like -f 120
##/home/cornejo/model/gridtest -lag 10 -x $i 10 -y $j 20 -f 127 -Name $k $arg
/home/cornejo/model/gridtest -lag 5 -x 20 10 -y 50 20 -f 127 -Name $k $arg
EOS

		done
	done
done

qstat -u cornejo
