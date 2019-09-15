echo "DISKSZ_MB_BCALM_GZ=$(echo $(du -k list_reads.unitigs.sset.gz | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')" >> global_stat 
gunzip list_reads.unitigs.sset.gz
echo "DISKSZ_MB_BCALM=$(echo $(du -k list_reads.unitigs.sset | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')" >> global_stat 

cat list_reads.unitigs.sset | awk -F=' ' '{print ">\n"$1}' > unitigs.mfc.txt
cat tipOutput.txt | tr "[" "a" | tr "]" "N" | awk -F=' ' '{print ">\n"$1}' > tip.mfc.txt
cat plainOutput.txt | awk -F=' ' '{print ">\n"$1}' > twoway.mfc.txt

$MFC=~/w/mfcompress/MFCompressC

$MFC unitigs.mfc.txt
DISKSZ_MB_BCALM_MFC=$(echo $(du -k unitigs.mfc.txt | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')
echo "DISKSZ_MB_BCALM_MFC=$DISKSZ_MB_BCALM_MFC" >> global_stat 

$MFC tip.mfc.txt
DISKSZ_MB_USTITCH_tip_MFC=$(echo $(du -k tip.mfc.txt.mfc | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')
echo "DISKSZ_MB_USTITCH_tip_MFC=$DISKSZ_MB_USTITCH_tip_MFC" >> global_stat 

$MFC twoway.mfc.txt
DISKSZ_MB_USTITCH_twoway_MFC=$(echo $(du -k twoway.mfc.txt.mfc | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')
echo "DISKSZ_MB_USTITCH_twoway_MFC=$DISKSZ_MB_USTITCH_twoway_MFC" >> global_stat 
