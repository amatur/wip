CC=g++
CFLAGS=-c -w -std=c++11 -O2

SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main.out

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o main.out incount.txt stitchedUnitigs.txt plainOutput.*

.SILENT:run

run:
	rm -f plainOutput.*
	#./main.out -i /Volumes/FAT32/data2019/phi11/list_reads.unitigs.fa -k 11 > myout.txt
	#./main.out -i /Volumes/exFAT/data2019/chol31/list_reads.unitigs.fa -f 1 -m 11 -k 31 > myout.txt
	./main.out -i /Volumes/exFAT/data2019/chol55/list_reads.unitigs.fa  -k 55 -f 1 -m 11 > old.txt

	cp stats.txt stats_sidedub.txt
	./main.out -i /Volumes/FAT32/data2019/chol55/list_reads.unitigs.fa  -k 55 -f 0 > myout.txt
	#./main.out -i /Volumes/FAT32/data2019/hum55/list_reads.unitigs.fa -k 55 
	#> myout.txt
	#./main.out -i /Users/Sherlock/cse566_2/exclude/bcalm-binaries-v2.2.1-Mac/bin/ecoli_genome21/list_reads.unitigs.fa -k 21 > myout.txt
	#./main.out -i /Users/Sherlock/cse566_2/exclude/bcalm-binaries-v2.2.1-Mac/bin/ecoli_genome102/list_reads.unitigs.fa -k 102 > myout.txt
	
	##PROFILE ONLY
	#./main.out -i /Volumes/FAT32/data2019/hum55/list_reads.unitigs.fa -k 55 -m 11 > myout.txt
	#./main.out -i /Users/Sherlock/cse566_2/exclude/bcalm-binaries-v2.2.1-Mac/bin/ecoli_genome21/list_reads.unitigs.fa -m 11 -k 21 > myout.txt
	#./main.out -i /Volumes/FAT32/data2019/chol31/list_reads.unitigs.fa -k 31 -m 11 > myout.txt
	
	cat stats_sidedub.txt
	cat stats.txt
	

	
	du -m plainOutput.txt | cut -f1 
	gzip plainOutput.txt
	du -m plainOutput.txt.gz | cut -f1 

	#./output-extract.sh

run45:
	./main.out -i /Volumes/FAT32/data2019/hum45/list_reads.unitigs.fa -k 45 

	