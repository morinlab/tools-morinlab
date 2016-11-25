# Get igv screenshot 

echo $@

# parse input parameters
set -- `getopt -n$0 -u -a --longoptions="track: region_file: archive: build: outfile: scriptdir: ftype: index: label: viewaspairs: view:" "h:" "$@"` || usage
[ $# -eq 0 ] && usage

while [ $# -gt 0 ]
do
    case "$1" in
                --track)	      	tracks+=",$2";shift;;  
		--region_file)		region_file=$2;shift;; 
		--build)		build=$2;shift;;
		--outfile)              outfile=$2;shift;;
		--scriptdir)            scriptdir=$2;shift;;
		--ftype)                ftypes+=",$2";shift;;
		--index)                indexfiles+=",$2";shift;;
		--label)                tracklabels+=",$2";shift;;
		--view)                 view+=",$2";shift;;
		--viewaspairs)          viewaspairs+=",$2";shift;;
                -h)       	 	shift;;
		--)       		shift;break;;
        -*)       		 	usage;;
        *)        		 	break;;            
    esac
    shift
done

# location of IGV installation
igvlocation="${scriptdir}/IGV_2.1.28"

# set preferences properly
echo -e "IGV.track.show.attribute.views=true
LAST_CHROMOSOME_VIEWED_KEY=All
DEFAULT_GENOME_KEY=${build}" > prefs.properties

#IGV.Bounds=0,0,$width,$height
#SHOW_SEQUENCE_TRANSLATION=true
#IGV.Frame.ExtendedState=0

# parse track, filetype and indexfile strings
tracks=${tracks/,/}     # remove leading comma from list
tracks=${tracks//,/ }   # replace commas with spaces
ftypes=${ftypes/,/}     # remove leading comma from list
indexfiles=${indexfiles/,/}       # remove leading comma from list
tracklabels=${tracklabels/,/}     # remove leading comma from list
tracklabels=${tracklabels// /_}   # replace spaces with underscores for track labels
view=${view/,/}     # remove leading comma from list
viewaspairs=${viewaspairs/,/}     # remove leading comma from list

# generate batch file recipe for igv
echo -e    "new 
genome ${build}" > recipe.txt

# add tracks to recipe
IFS=',' read -a ftypearray <<< "$ftypes"
IFS=',' read -a indexarray <<< "$indexfiles"
IFS=',' read -a labelarray <<< "$tracklabels"
IFS=',' read -a viewarray <<< "$view"
IFS=',' read -a viewaspairsarray <<< "$viewaspairs"
count=0
for track in $tracks
do 

    # make symlink to file with the proper extension (or IGV will not work)
    ftype=${ftypearray[count]}
    label=${labelarray[count]}
    ltrack="${label}.${ftype}"
    ln -s ${track} ${ltrack}

    # link the index file too if needed
    if [[ $ftype == "bam" ]]
    then
        indexltrack="${label}.bam.bai"
        ln -s ${indexarray[count]} ${indexltrack}
    fi

    # add track to recipe
    echo -e "load ${ltrack}" >>recipe.txt

    # set view for track (collapsed/expanded/squished)
    echo -e "${viewarray[count]} ${ltrack}" >> recipe.txt
    
    # set to view as pairs if requested
    if [[ ${viewaspairsarray[count]} == "yes" ]]
    then
        echo -e "viewaspairs ${ltrack}" >> recipe.txt
    fi
    
    
    count=$[$count+1]

done

#use awk to convert bed-style format into required chr:s-e format

echo -e "snapshotDirectory outputs
maxPanelHeight 3000" >> recipe.txt
awk '{$1;$2;$3;$out = $1":"$2-40"-"$3+40;print $out}' < $region_file > regions_reformat.txt
while IFS= read -r line
do
    echo -e "goto $line
    sort base
    snapshot $line.png" >> recipe.txt
done < "regions_reformat.txt"
echo "
exit" >> recipe.txt

# for debugging purposes:
echo -e "\nIGV batch recipe is:\n"
cat recipe.txt

# setup X environment
echo -e "\nsetting up Xvfb:\n"
Xvfb :12 2>&1 &
export DISPLAY=:12

# run igv
echo -e "\nrunning IGV:\n"
igv.sh -b recipe.txt -o prefs.properties 2>&1

