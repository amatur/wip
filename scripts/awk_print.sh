printf "" > mystat.txt
cat global_stat | awk -F= '$1=="DISKSZ_MB_USTITCH_tip"{print $2}' >> mystat.txt
cat global_stat | awk -F= '$1=="DISKSZ_MB_USTITCH_tip_GZ"{print $2}'  >> mystat.txt
cat global_stat | awk -F= '$1=="DISKSZ_MB_USTITCH_twoway"{print $2}'  >> mystat.txt
cat global_stat | awk -F= '$1=="DISKSZ_MB_USTITCH_twoway_GZ"{print $2}'  >> mystat.txt

awk 'BEGIN{ORS="\t"}1' mystat.txt 
