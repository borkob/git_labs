#/bin/bash
targetReached=0;
coordInit="11101011100"
maxWalkLength=0
minWalkLength=10000000000
echo ${coordInit}
rm -f walkLength.txt
rm -f walks32.txt
RANDOM=777
program=./lssRRts
part="3"
if [ "$part" = "disabled" ] || [ "$part" = "1" ]
then
	for i in {1..32} 
	do
		walkLength=0
		seed=$RANDOM
		$program 21 ${seed} 300 1 26 -walk -coordInit ${coordInit} | awk '{if($2 == 0) print $0}' > walk.txt
		echo [${i}/1] $program 21 ${seed} 300 1 26 -walk -coordInit ${coordInit} >>  walks32.txt 
		cat walk.txt >> walks32.txt
		walkLength=`cat walk.txt | wc -l`
		walkLength=$(( $walkLength - 1 ))
		echo [${i}/1] $program 21 ${seed} 300 1 26 -walk -coordInit ${coordInit}
		tr=`cat walk.txt | grep targetReached | head -n 1`
	        if [ -z "$tr" ]
        	then
                	tr=`cat walk.txt | tail -n 1 | awk '{print $4}'`
			if [ "$tr" -eq "26" ]
	                then
        	                walkLength=0
                	else
                        	continue
	                fi
        	else
                	walkLength=`echo $tr |  awk '{print $1}'`
	        fi
		walkLength=`echo $tr |  awk '{print $1}'`
		echo [${i}/1] "Walk length:" $walkLength
		echo ${walkLength} >> walkLength.txt
		targetReached=$(( $targetReached + 1 ));
		echo [${i}/1] "Target reached:" $targetReached
		if [ $walkLength -lt $minWalkLength ]
		then
			minWalkLength=${walkLength}
		fi
		if [ $walkLength -gt $maxWalkLength ]
		then
			maxWalkLength=${walkLength}
		fi
	done;
	echo Walk length [${minWalkLength},${maxWalkLength}]
	echo Target reached: $targetReached 
	R --slave --no-save < results.R
fi

if [ "$part" = "disabled" ] || [ "$part" = "2" ]
then
	rm -f walks1000.txt walks1000.txt
	targetReached=0;
	rm -f walkLength.txt
	for i in {1..1000} 
	do
		walkLength=0
		seed=$RANDOM
		$program 21 ${seed} 300 1 26 -walk | awk '{if($2 == 0) print $0}' > walk.txt
		echo [${i}/2] $program 21 ${seed} 300 1 26 -walk >> walks1000.txt
		cat walk.txt >> walks1000.txt
		walkLength=`cat walk.txt | wc -l`
		walkLength=$(( $walkLength - 1 ))
		echo [${i}/2] $program 21 ${seed} 300 1 26 -walk 
		tr=`cat walk.txt | grep targetReached | head -n 1`
		if [ -z "$tr" ]
		then
			tr=`cat walk.txt | tail -n 1 | awk '{print $4}'`
			if [ "$tr" -eq "26" ]
			then
				walkLength=0
			else
				continue
			fi
		else
			walkLength=`echo $tr |  awk '{print $1}'`
		fi
		echo [${i}/2] "Walk length:" $walkLength
		echo $walkLength >> walkLength.txt
		targetReached=$(( $targetReached + 1 ));
		echo [${i}/2] "Target reached:" $targetReached
		if [ $walkLength -lt $minWalkLength ]
		then
			minWalkLength=$walkLength
		fi
		if [ $walkLength -gt $maxWalkLength ]
		then
			maxWalkLength=$walkLength
		fi
	done;
	echo $program
	echo Walk length [${minWalkLength},${maxWalkLength}]
	echo Target reached: $targetReached 
	R --slave --no-save < results.R
fi

L=27
T=37
if [ "$part" = "all" ] || [ "$part" = "3" ]
then
	targetReached=0;
	rm -f walkLength.txt walks2_1000.txt
	for i in {1..10000} 
	do
		walkLength=0
		seed=$RANDOM
		$program $L ${seed} 300 1 $T -walk > walk.txt
		#echo [${i}/3] $program $L ${seed} 300 1 $T -walk >> walks2_1000.txt
		#cat walk.txt >> walks2_1000.txt
		walkLength=`cat walk.txt | wc -l`
		walkLength=$(( $walkLength - 1 ))
		echo [${i}/3] $program $L ${seed} 300 1 $T -walk 
		tr=`cat walk.txt | grep targetReached | head -n 1`
		if [ -z "$tr" ]
		then
                	tr=`cat walk.txt | tail -n 1 | awk '{print $4}'`
        	        if [ "$tr" -eq "$T" ]
                	then
                        	walkLength=0
	                else
				echo "Error: target is not reached !" 
				exit 1
                	fi
		else
			walkLength=`echo $tr |  awk '{print $1}'`
		fi
		echo [${i}/3] "Walk length:" $walkLength
		echo $walkLength >> walkLength.txt
		targetReached=$(( $targetReached + 1 ));
		echo [${i}/3] "Target reached:" $targetReached
		if [ $walkLength -lt $minWalkLength ]
		then
			minWalkLength=$walkLength
		fi
		if [ $walkLength -gt $maxWalkLength ]
		then
			maxWalkLength=$walkLength
		fi
	done;
	echo $program
	echo Walk length [${minWalkLength},${maxWalkLength}]
	echo Target reached: $targetReached 
	R --slave --no-save < results.R
fi
program=./lssMAts
if [ "$part" = "all" ] || [ "$part" = "4" ]
then
	targetReached=0;
	rm -f walkLength.txt walks2_1000.txt
	for i in {1..10000} 
	do
		walkLength=0
		seed=$RANDOM
		$program 21 ${seed} 300 1 26 -walk > walk.txt
		#echo [${i}/3] $program 21 ${seed} 300 1 26 -walk >> walks2_1000.txt
		#cat walk.txt >> walks2_1000.txt
		walkLength=`cat walk.txt | wc -l`
		walkLength=$(( $walkLength - 1 ))
		echo [${i}/3] $program 21 ${seed} 300 1 26 -walk 
		tr=`cat walk.txt | grep targetReached | head -n 1`
		if [ -z "$tr" ]
		then
                	tr=`cat walk.txt | tail -n 1 | awk '{print $4}'`
        	        if [ "$tr" -eq "26" ]
                	then
                        	walkLength=0
	                else
				echo "Error: target is not reached !" 
				exit 1
                	fi
		else
			walkLength=`echo $tr |  awk '{print $1}'`
		fi
		echo [${i}/3] "Walk length:" $walkLength
		echo $walkLength >> walkLength.txt
		targetReached=$(( $targetReached + 1 ));
		echo [${i}/3] "Target reached:" $targetReached
		if [ $walkLength -lt $minWalkLength ]
		then
			minWalkLength=$walkLength
		fi
		if [ $walkLength -gt $maxWalkLength ]
		then
			maxWalkLength=$walkLength
		fi
	done;
	echo $program
	echo Walk length [${minWalkLength},${maxWalkLength}]
	echo Target reached: $targetReached 
	R --slave --no-save < results.R
fi
