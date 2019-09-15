
#include <cmath>
#include<cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdint.h>
#include <unistd.h>

int K = 11;
using namespace std;
int get_data(int K, string UNITIG_FILE= "/Users/Sherlock/Documents/bcl/bcl/tipOutput.txt"){
    ifstream unitigFile;
    unitigFile.open(UNITIG_FILE);
    
    ofstream outFile;
    outFile.open("tip.fa");
    
    string line;
    bool startPrefCut = false;
    bool startSufCut = false;
    string pref;
    string suf;
    int walkid = 0;
    
    while (getline(unitigFile, line)) {
        cout<<line<<endl;
        if (line.empty() || line.substr(0, 1).compare(">") == 0) {
            
        } else {
            walkid++;

            startPrefCut = false;
            startSufCut = false;
            
            string tip = "";
            string sbuf = "";
            string lastk = "";
            string pref = "";
            string suf = "";
            
            cout<<walkid-1<<" "<<line[0]<<endl;
            for(int i = 0; i<line.length(); i++){
                if(line[i]=='A'|| line[i]=='C' || line[i]=='G' || line[i]=='T'){
                    if(startPrefCut){
                        tip += line[i];
                    }else if(startSufCut){
                        tip += line[i];
                    }else{
                        sbuf +=   line[i];
                        if(lastk.length() < K - 1){
                            lastk += line[i];
                        }else{
                            lastk = lastk.substr(1, K-2) + line[i];
                        }
                    }
                }else if(line[i]==']'){
                    if(startPrefCut){ //already has one
                        startPrefCut = false;
                        outFile<< "> pref: " << " "  << endl;
                        outFile<< (pref + tip) << endl;
                        tip = "";
                    }else if(!startPrefCut){ //prefix is cut starts
                        startPrefCut = true;
                        pref = lastk;
                    }
                }else if(line[i]=='['){ //suffix is cut
                    if(startSufCut){ //already has one
                        startSufCut = false;
                        outFile<< ">suf: " << " " << endl;
                        outFile<< (tip + suf) << endl;
                        tip = "";
                    }else if(!startSufCut){
                        startSufCut = true;
                        suf = lastk;
                    }
                }
            }
            //print the tips before printing the contig
            outFile<<">"<<walkid-1<<" : "<<"\n";
            outFile<<sbuf << endl;
        }
    }
    
    unitigFile.close();
    outFile.close();

    cout << "Complete conversion." << endl;
    return EXIT_SUCCESS;
}

int main(int argc, char** argv){
    const char* stringValue = "" ;
    int c ;
    
    ///*
    while( ( c = getopt (argc, argv, "i:k:m:d:f:") ) != -1 )
    {
        switch(c)
        {
            case 'i':
                if(optarg) stringValue = optarg;
                break;
            case 'k':
                if(optarg) {
                    K = std::atoi(optarg) ;
                    if(K<=0){
                        fprintf(stderr, "Error: Specify a positive k value.\n");
                        exit(EXIT_FAILURE);
                    }
                }else{
                    fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                            argv[0]);
                    exit(EXIT_FAILURE);
                }
                break;
            default: //
                fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                        argv[0]);
                exit(EXIT_FAILURE);
                
        }
    }
    
    if(K==0 || strcmp(stringValue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
   //    */
    
    get_data(K);
    return 0;
}
