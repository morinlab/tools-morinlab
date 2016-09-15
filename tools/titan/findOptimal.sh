val=1;                                                                                                                                                                
index=0;
for i in $(eval echo "{$1..$2}");
    do newval=$( grep ".*S_Dbw validity index (Both).*" ./parameters/samp${i}.txt | cut -f2);
    if [ $(echo "$val > $newval" | bc) -eq 1 ];
        then val=$newval; index=$i;
    fi;
done;
echo -e "$val\t$index"
