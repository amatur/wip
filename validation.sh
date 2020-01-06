K=$2

cd /Users/Sherlock/Documents/bcl/bcl/

/Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk  -file $1  -kmer-size $K -abundance-min 1 
#~ /Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk  -file /Users/Sherlock/Documents/bcl/bcl/stitchedUnitigs.fa  -kmer-size $K -abundance-min 1 
#~ /Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk  -file /Users/Sherlock/Documents/bcl/bcl/stitchedUnitigsFwd.fa  -kmer-size $K -abundance-min 1 
/Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk  -file /Users/Sherlock/Documents/bcl/bcl/stitchedUnitigsTip.fa  -kmer-size $K -abundance-min 1 


/Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk2ascii -file list_reads.unitigs.h5 -out output-bcalm.txt
#~ /Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk2ascii -file stitchedUnitigs.h5 -out output-my.txt
#~ /Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk2ascii -file stitchedUnitigsFwd.h5 -out output-fwd.txt
/Volumes/exFAT/work/dsk-v2.3.0-bin-Darwin/bin/dsk2ascii -file stitchedUnitigsTip.h5 -out output-tip.txt


cat output-bcalm.txt | awk '{ sum += $2 } END { printf "BCALM2 = %d\n", sum }'
#~ cat output-my.txt | awk '{ sum += $2 } END { printf "USTITCH = %d\n", sum }'
#~ cat output-fwd.txt | awk '{ sum += $2 } END { printf "USTITCH FWD = %d\n", sum }'
cat output-tip.txt | awk '{ sum += $2 } END { printf "USTITCH TIP = %d\n", sum }'

#sort output-my.txt > a.txt
#sort output-bcalm.txt > b.txt

#diff a.txt b.txt
