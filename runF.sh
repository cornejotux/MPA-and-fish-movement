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
			walltime="40:00:00"
			arg="-Int -SST -Cur"
			echo "1 $arg"
		else
			walltime="30:00:00"
			arg="-Int -SST -NoCur"
			echo "2 $arg"
		fi
	else
		walltime="20:00:00"
		arg="-Int -NoSST -NoCur"	
		echo "3 $arg"
	fi
else
	walltime="10:00:00"
	arg="-NoSST -NoCur -NoInt"
	echo "4 $arg"
fi
## This are the annual fishing rates x 100, because the shell do not work with decimals
scenarios=("50" "55" "60" "65" "70" "75" "80" "85" "90" "95" "100" "105" "110" "115" "120" "125" "130")	

for i in ${scenarios[@]}
do
    cat <<EOS | qsub -
#!/bin/sh

#PBS -N $1_$2_$3_$i
#PBS -l nodes=1,walltime=$walltime

cd /home/cornejo/model/
/home/cornejo/model/gridtest -MPA 15 $arg -F $i
EOS
done

qstat -u cornejo
