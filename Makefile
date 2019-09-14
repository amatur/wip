CC=g++
CFLAGS=-c -w -std=c++11 -O2

SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main.out
BCALMFILE=/Volumes/exFAT/data2019/phi/11/list_reads.unitigs.fa
K=11

all: $(SOURCES) $(EXECUTABLE) decoderd

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

decoderd: decoder.cpp
	$(CC) $(CFLAGS) decoder.cpp -o decoder.o
	$(CC) decoder.o -o decoder.out
	rm decoder.o

clean:
	rm -f *.o main.out incount.txt stitchedUnitigs.txt plainOutput.* stitched* plain* *.fa *.txt

.SILENT:run

run:
	#rm -f *.txt *.fa
	./main.out -i  $(BCALMFILE) -k $(K) -f 0 -m 15 > myout.txt
	./decoder.out -i tipOutput.txt -k $(K) > decot.txt

	open decot.txt
	/Volumes/exFAT/work/validation.sh $(BCALMFILE) $(K) 
	#./main.out -i /Volumes/exFAT/data2019/chol31/list_reads.unitigs.fa -f 1 -m 11 -k 31 > myout.txt
	#./main.out -i /Volumes/exFAT/data2019/chol/31/list_reads.unitigs.fa  -k 31 -f 1 -m 0 > myout0.txt
	#./main.out -i /Volumes/exFAT/data2019/chol/31/list_reads.unitigs.fa  -k 31 -f 1 -m 10 > myout10.txt
	#./main.out -i /Volumes/exFAT/data2019/staph31/list_reads.unitigs.fa  -k 31 -f 1 -m 0 > myout15.txt
	#cat myout.txt 

	#cp stats.txt stats_sidedub.txt

	#cat stats_sidedub.txt
	#cat stats.txt
	
	#open myout0.txt
	#open myout10.txt
	open myout.txt

	#./validation.sh /Volumes/exFAT/data2019/chol/55/list_reads.unitigs.fa 55  
	
	#du -m plainOutput.txt | cut -f1 
	#gzip plainOutput.txt
	#du -m plainOutput.txt.gz | cut -f1 
