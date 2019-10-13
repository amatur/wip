import sys
filepath = (sys.argv[1])
kmer_size = int(sys.argv[2])
fout = open("kmers_serially.fa", "w+")

with open(filepath) as fp:
    cnt = 0
    for string in fp:
        if(string[0]=='>'):
            continue
        seqlen = len(string.strip())
        for i in range(0, seqlen - kmer_size + 1):
            string = list(string)
            kmer = ''.join(string[i:i + kmer_size])
            fout.write(">\n%s\n" % kmer)
fout.close()

