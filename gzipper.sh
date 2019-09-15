DISKSZ_MB_USTITCH_tip=$(echo $(du -k tipOutput.txt | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')
echo "DISKSZ_MB_USTITCH_tip=$DISKSZ_MB_USTITCH_tip" >> global_stat 
rm -f tipOutput.txt.gz 
gzip tipOutput.txt
DISKSZ_MB_USTITCH_tip_GZ=$(echo $(du -k tipOutput.txt.gz | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')
echo "DISKSZ_MB_USTITCH_tip_GZ=$DISKSZ_MB_USTITCH_tip_GZ" >> global_stat 


echo "DISKSZ_MB_USTITCH_twoway=$(echo $(du -k plainOutput.fa | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')" >> global_stat 
rm -f plainOutput.fa.gz 
gzip plainOutput.fa
echo "DISKSZ_MB_USTITCH_twoway_GZ=$(echo $(du -k plainOutput.fa.gz | cut -f 1) | awk '{printf "%.2f\n",$1/1024.0}')" >> global_stat 
