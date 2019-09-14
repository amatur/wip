
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

int K = 11;

using namespace std;



int get_data(int K, string UNITIG_FILE= "/Users/Sherlock/Documents/bcl/bcl/tipOutput.txt"){
    ifstream unitigFile;
    unitigFile.open(UNITIG_FILE);
    //    if(!unitigFile.good()){
    //        fprintf(stderr, "Error: File named \"%s\" cannot be opened.\n", UNITIG_FILE.c_str());
    //        exit(EXIT_FAILURE);
    //    }
    
    ofstream outFile;
    outFile.open("tip.fa");
    
    string line;
    
    getline(unitigFile, line); // a ">"
    
    int it = 0;
    bool startLeft = false;
    bool startRight = false;
    //bool endLeft = false;
    //bool endRight = false;
    char cutType = 'n';
    string pref;
    string suf;
    int tipIt = 0;
    int tipStartIndex = -1;
    
    //char tipbuffer[2000000];
    //char sbuffer[2000000];
    
    int walkid = 0;
    
    while (getline(unitigFile, line)) {
        cout<<line<<endl;
        if (line.empty() || line.substr(0, 1).compare(">") == 0) {
            //unitig_struct.sequence = unitig_struct.sequence + line;
            //unitigs.push_back(unitig_struct);
        } else {
            walkid++;
            if(walkid==32){
                
            }
            
           startLeft = false;
           startRight = false;
            
            string tip = "";
            string sbuf = "";
            tipIt = 0;
            it = 0;
            string lastk = "";
            
            cout<<walkid-1<<" "<<line[0]<<endl;
            for(int i = 0; i<line.length(); i++){
                if(line[i]=='A'|| line[i]=='C' || line[i]=='G' || line[i]=='T'){
                    if(startLeft){
                        tip += line[i];
                    }else if(startRight){
                        tip += line[i];
                    }else{
                        //sbuffer[it++] = line[i];
                        sbuf +=   line[i];
                        if(lastk.length() < K - 1){
                            lastk += line[i];
                        }else{
                            lastk = lastk.substr(1, K-2) + line[i];
                        }
                        
                    }
                }else if(line[i]==']'){
                    if(i<K)cout<<i<<endl;
                    if(startLeft){ //already has one
                        startLeft = false;
                        //tipbuffer[tipIt] = '\0';
                        
                        //string tt(tipbuffer);
                        string tt = tip;
                        outFile<< "> pref: " <<tipIt<< " "  << endl;
                        tipIt = 0;
                        outFile<< (pref + tt);
                        outFile<<"\n";
                        
                    }
                    if(!startLeft){
                        //prefix is cut
                        tipIt = 0;
                        startLeft = true;
                        cutType = 'p';
                        //if(line[i-1]=='A'|| line[i-1]=='C' || line[i-1]=='G' || line[i-1]=='T'){
                        //    tipStartIndex = i;
                        //}
                        //int c = tipStartIndex - (K - 1);
                        // pref = line.substr(tipStartIndex - (K - 1), (K-1));
                        pref = lastk;
                        //                        pref = "";
                        //                        for(int xx = 0; xx < (K-1); xx++){
                        //                            pref += sbuffer[it-K+1+xx];
                        //                        }
                        // cout<<pref<<endl;
                    }
                }else if(line[i]=='['){ //suffix is cut
                    if(i<31)cout<<i<<endl;
                    if(startRight){ //already has one
                        //endRight = true;
                        startRight = false;
                        //tipbuffer[tipIt] = '\0';
                        
                        //string tt(tipbuffer);
                        string tt = tip;
                        outFile<< ">suf " <<tipIt<< " " << endl;
                        
                        tipIt = 0;
                        outFile<< (tt + suf);
                        outFile<<"\n";
                    }
                    if(!startRight){
                        tipIt = 0;
                        startRight = true;
                        cutType = 's';
                        //if(line[i-1]=='A'|| line[i-1]=='C' || line[i-1]=='G' || line[i-1]=='T'){
                        //    tipStartIndex = i;
                        //}
                        suf = lastk;
                        //
                        //                        suf = "";
                        //                        for(int xx = 0; xx < (K-1); xx++){
                        //                            suf += sbuffer[it-K+1+xx];
                        //                        }
                        //cout<<suf<<endl;
                        
                    }
                }
            }
            //print the tips before printing the contig
            //sbuffer[it]='\0';

                //string oo(sbuffer);
                outFile<<">"<<walkid-1<<" : "<<it<<"\n";
                outFile<<sbuf;
                outFile<<"\n";
            
            
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
    //system("cp tipOutput.txt tip2.txt");
    //system("gzip tipgz.txt");
    
    
    
    
    return 0;
}
