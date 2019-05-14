if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <readfile> <proportion_reads_long>"
    exit 1
fi

readfile=$1
prop=$2

filetype='fastq'
firstchar=`head -c 1 $readfile`
if  [ "$firstchar" = ">" ]
then
    filetype='fasta'
fi

if [ "$filetype" = "fastq" ]
then
    cat $readfile | awk '{if(NR%4==2) print length($1)}' | sort -n -r | uniq -c > read_length.txt
else
    cat $readfile | awk '{if(NR%2==0) print length($1)}' | sort -n -r | uniq -c > read_length.txt
fi
numreads=`awk '{s+=$1} END {print s}' read_length.txt`
echo 'NumReads: '$numreads
keptreads=`echo "$(echo "scale=4; $prop*$numreads" | bc)" | awk '{printf("%d\n",$1)}'`
echo 'KeptReads: '$keptreads

length_cutoff=`awk -v kr="$keptreads" 'BEGIN {s = 0; last = 0;} \
        { \
            s += $1; \
            if (s >= kr) { \
                print last; \
                exit; \
            } \
            last = $2; \
        }' read_length.txt`

echo 'LengthCutoff: '$length_cutoff

srf=$readfile'.short'
srf=''
lrf=$readfile'.long'

echo 'LongReadFile: '$lrf
echo 'ShortReadFile: '$srf

awk -v lc="$length_cutoff" -v srf="$srf" -v lrf="$lrf" -v ft="$filetype" 'BEGIN {FS = "\t" ; OFS = "\n"} 
        { \
            header = $0 ; \
            getline seq ; \
            if (ft == "fastq") { \
                getline qheader ; \
                getline qseq ; \
                if (length(seq) >= lc) { \
                    print header, seq, qheader, qseq > lrf; \
                } \
                else if(length(srf) > 0) { \
                   print header, seq, qheader, qseq > srf; \
                } \
            } \
            else { \
               if (length(seq) >= lc) { \
                    print header, seq > lrf; \
                } \
                else if(length(srf) > 0) { \
                   print header, seq > srf; \
                } \
            } \
        }' < $readfile
