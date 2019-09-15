import os
import sys
import subprocess

if(len(sys.argv)<2):
    os.chdir("/Users/Sherlock/Documents/bcl/bcl")
else:
    os.chdir(sys.argv[1])

statDict = {}
filepath = 'global_stat'
with open(filepath) as fp:
    for line in fp:
        x = (line.strip().split("="))
        statDict[x[0]] = x[1]
# print(statDict)
fields = "K	N_KMER	C_BCALM	C_LB	C_USTITCH_twoway	C_USTITCH_tip	V_BCALM	V_LB	V_USTITCH_twoway	V_USTITCH_tip	DISKSZ_MB_BCALM	DISKSZ_MB_USTITCH_twoway	DISKSZ_MB_USTITCH_tip	DISKSZ_MB_BCALM_GZ	DISKSZ_MB_USTITCH_twoway_GZ	DISKSZ_MB_USTITCH_tip_GZ".split()
print(fields)
for key, value in statDict.items() :
    print(key, " ", end='\t')


os.remove("sub_stat.txt")
fout = open("sub_stat.txt", "w+")
for i in fields:
    if(i in statDict):
        fout.write("%s\t" % (statDict[i]))
    else:
        fout.write("%s\t" % 0)
    # print("%s\t" % (dict[i]))
fout.close()
subprocess.call(['open', "sub_stat.txt"])
