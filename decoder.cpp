
#include <cmath>
#include<cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdint.h>
#include <unordered_set>
#include <map>
#include <set>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <list>
#include <stack>
#include <unordered_map>
#include <utility>
#include <queue>
#include <deque>
#include <unistd.h>
#include <tuple>

int K = 31;

using namespace std;



int get_data(int K){
    ifstream unitigFile;
    unitigFile.open("/Users/Sherlock/Documents/bcl/bcl/decoder/decoder/tipOutput.txt");
    
    ofstream outFile;
    outFile.open("/Users/Sherlock/Documents/bcl/bcl/decoder/decoder/tipVali.fa");
    
    string line;
    
    getline(unitigFile, line); // a ">"
    
    size_t it = 0;
    bool startLeft = false;
    bool startRight = false;
    bool endLeft = false;
    bool endRight = false;
    char cutType = 'n';
    string pref;
    string suf;
    int tipIt = 0;
    int tipStartIndex = -1;
    
    char tipbuffer[2000000];
    char sbuffer[2000000];
    
    while (getline(unitigFile, line)) {
        //cout<<line;
        if (line.substr(0, 1).compare(">") == 1) {
            //unitig_struct.sequence = unitig_struct.sequence + line;
            //unitigs.push_back(unitig_struct);
        } else {
            it = 0;
            size_t slen = line.length();
            for(size_t i = 0; i<slen; i++){
                if(line[i]=='A'|| line[i]=='C' || line[i]=='G' || line[i]=='T'){
                    if(startLeft && !endLeft){
                        tipbuffer[tipIt++] = line[i];
                    }else if(!startLeft && endLeft){
                        tipbuffer[tipIt] = '\0';
                        string tt(tipbuffer);
                        outFile<< ">\n";
                        outFile<< (pref + tt);
                        outFile<<"\n";
                        endLeft = false;
                    }else if(startRight && !endRight){
                        tipbuffer[tipIt++] = line[i];
                    }else if(!startRight && endRight){
                        tipbuffer[tipIt] = '\0';
                        string tt(tipbuffer);
                        outFile<< ">\n";
                        outFile<< (tt + suf);
                        outFile<<"\n";
                        endRight = false;
                    }else{
                        sbuffer[it++] = line[i];
                    }
                }else if(line[i]==']'){
                    if(startLeft && !endLeft){ //already has one
                        endLeft = true;
                        startLeft = false;
                    }
                    if(!startLeft){
                        //prefix is cut
                        tipIt = 0;
                        startLeft = true;
                        cutType = 'p';
                        tipStartIndex = i;
                        pref = line.substr(tipStartIndex - K - 1, (K-1));
                    }
                }else if(line[i]=='['){ //suffix is cut
                    if(startRight && !endRight){ //already has one
                        endRight = true;
                        startRight = false;
                    }
                    if(!startRight){
                        tipIt = 0;
                        startRight = true;
                        cutType = 's';
                        tipStartIndex = i;
                        suf = line.substr(tipStartIndex - K - 1, (K-1));
                    }
                }
            }
            //print the tips before printing the contig
            sbuffer[it]='\0';
            //cout<<it;
            
            string oo(sbuffer);
            outFile<<">\n";
            outFile<<oo;
            outFile<<"\n";
        }
    }
    
    
    unitigFile.close();
    outFile.close();
    
    cout << "Complete reading input." << endl;
    
    return EXIT_SUCCESS;
}

int main(){
    get_data(K);
    //system("cp tipOutput.txt tip2.txt");
    //system("gzip tipgz.txt");
    return 0;
}
