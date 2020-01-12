//
//  global.h
//  bcl
//
//  Created by Amatur Rahman on 24/11/19.
//  Copyright © 2019 psu. All rights reserved.
//

#ifndef global_h
#define global_h

bool DEBUGMODE = false;

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>
#include <stdint.h>
#include <unordered_set>
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
#include <set>
#include <map>

using namespace std;

//int K = 11;
//string UNITIG_FILE = "/Volumes/exFAT/data2019/staphsub/11/list_reads.unitigs.fa";

//int K = 31;
//string UNITIG_FILE = "/Volumes/exFAT/data2019/chol/31/list_reads.unitigs.fa";

int K = 11;
string UNITIG_FILE = "/Users/Sherlock/amaturWS/data2019/staphsub/"+to_string(K)+"/list_reads.unitigs.fa";
//string UNITIG_FILE = "/Users/Sherlock/amaturWS/data2019/staphsub/11/list_reads.unitigs.fa";

enum DEBUGFLAG_T { NONE = 0,  UKDEBUG = 0, VERIFYINPUT = 1, INDEGREEPRINT = 2, DFSDEBUGG = 3, PARTICULAR = 4, NODENUMBER_DBG = 5, OLDNEWMAP = 9, PRINTER = 10, SINKSOURCE = 12};

enum ALGOMODE_T { BASIC = 0, INDEGREE_DFS = 1, INDEGREE_DFS_1 = 2, OUTDEGREE_DFS = 3, OUTDEGREE_DFS_1 = 4, INDEGREE_DFS_INVERTED = 5, PLUS_INDEGREE_DFS = 6, RANDOM_DFS = 7, NODEASSIGN = 8, SOURCEFIRST = 9, TWOWAYEXT = 10, PROFILE_ONLY = 11, EPPRIOR=12, GRAPHPRINT = 13, TIGHTUB = 14, BRACKETCOMP = 15, ONEWAYABSORPTION=16, ONEWAYABSORPTION_UNTESTED = 17};

bool FLG_NEWUB = true;

DEBUGFLAG_T DBGFLAG = NONE; //NODENUMBER_DBG
ALGOMODE_T ALGOMODE = ONEWAYABSORPTION;

bool MODE_WALK_UNION = (ALGOMODE == TWOWAYEXT);
bool MODE_ABSORPTION_TIP = (ALGOMODE == BRACKETCOMP);
bool MODE_ABSORPTION_NOTIP = (ALGOMODE == ONEWAYABSORPTION || ALGOMODE == ONEWAYABSORPTION_UNTESTED);


string mapmode[] = {"basic", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "twoway", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "tip", "one_way_absorption", "one_way_absorption_not_tested"
};
string modefilename[] = {"Fwd", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "Tip", "OneAbsorption", "one_way_absorption_not_tested"
};

/******
 FILENAMES
 
 */

namespace MyTypes
{
    typedef tuple<int,int,int, int> fourtuple; // uid, walkid, pos, isTip
    // etc.
}

typedef struct {
    int serial = -1;
    int startPosWithKOverlap;
    int endPosWithKOVerlap;
    bool isWalkEnd = false;
    int pos_in_walk = -100;
    int finalWalkId = -1; // renders some walkId as invalid
    int isTip = 0;
} new_node_info_t;

typedef struct {
    int serial;
    string sequence;
    int ln;
    int kc;
    float km;
} unitig_struct_t;
vector<unitig_struct_t> unitigs;


typedef struct {
    //1 means +, 0 means -
    bool left;
    bool right;
    int toNode;
} edge_t;

typedef struct {
    edge_t edge;
    int fromNode;
} edge_both_t;

typedef struct {
    edge_t edge;
    int kmerStartIndex;
    int kmerEndIndex;
} newEdge_t;

bool sort_by_walkId (const MyTypes::fourtuple &lhs, const MyTypes::fourtuple &rhs){
    return get<1>(lhs) < get<1>(rhs);
}
bool sort_by_pos (const MyTypes::fourtuple &lhs, const MyTypes::fourtuple &rhs){
    return get<2>(lhs) < get<2>(rhs);
}
bool sort_by_tipstatus (const MyTypes::fourtuple &lhs, const MyTypes::fourtuple &rhs){
    return get<3>(lhs) < get<3>(rhs);
}

struct node_sorter {
    int node;
    int sortkey;
    //bool operator() (struct node_sorter  i, struct node_sorter  j) { return (i.sortkey<j.sortkey);}
};
bool sort_by_key (struct node_sorter i, struct node_sorter j) { return (i.sortkey<j.sortkey); }
bool sort_by_key_inverted (struct node_sorter i, struct node_sorter j) { return (i.sortkey>j.sortkey); }

int* global_indegree;
int* global_outdegree;
int* global_plusindegree;
int* global_plusoutdegree;
int* global_issinksource;
int* global_priority;
new_node_info_t* oldToNew;
bool* nodeSign;
//absorber-global
bool* obsoleteWalkId;
int countNewNode = 0; // number of walks by ust-onewaymerge

map<pair <int, int>, int> inOutCombo;
vector<vector<edge_t> > adjList;
vector<vector<edge_t> > reverseAdjList;

map<int, string> newSequences;
map<int, string> newNewSequences; //int is the unitig id (old id)
set<int> newNewMarker;

vector<list<int> > newToOld;
vector<int> walkFirstNode; //given a walk id, what's the first node of that walk
unordered_map<int, vector<edge_t> > sinkSrcEdges; //int is the unitig id (old id)


//connected component count graph
vector<set<int> > ccAdjList;
int absorbGraphNumCC=0;


int C_ustitch = 0;
int C_twoway_ustitch = 0;
int C_tip_ustitch = 0;
int C_oneabsorb = 0;
int C_oneabsorb_ACGT = 0;
int C_oneabsorb_brackets = 0;
int C_oneabsorb_plusminus = 0;


int V_ustitch = 0;
int V_twoway_ustitch = 0;
int V_tip_ustitch = 0;
int V_oneabsorb = 0;

int isolated_node_count = 0;
int sink_count = 0;
int source_count = 0;
int sharedparent_count = 0;
int sharparentCntRefined = 0;
int onecount = 0;

inline string plus_strings(const string& a, const string& b, size_t kmersize) {
    if (a == "") return b;
    if (b == "") return a;
    string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
    return ret;
}

inline string pref(const string& b, size_t kmersize) {
    return b.substr(0, (kmersize - 1));
}
inline string cutPref(const string& b, size_t kmersize) {
    return b.substr(kmersize - 1, b.length() - (kmersize - 1));
}
inline string suf(const string& b, size_t kmersize) {
    return b.substr(b.length() - (kmersize - 1), (kmersize - 1) );
}
inline string cutSuf(const string& b, size_t kmersize) {
    return b.substr(0, b.length() - (kmersize - 1));
}


inline string delSpaces(string &str) {
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

inline bool charToBool(char c) {
    if (c == '+') {
        return true;
    } else {
        if (c != '-') cout << "ERRRRRROOR!" << endl;
        return false;
    }
}

inline string reverseComplement(string base) {
    size_t len = base.length();
    char* out = new char[len + 1];
    out[len] = '\0';
    for (int i = 0; i < len; i++) {
        if (base[i] == 'A') out[len - i - 1] = 'T';
        else if (base[i] == 'C') out[len - i - 1] = 'G';
        else if (base[i] == 'G') out[len - i - 1] = 'C';
        else if (base[i] == 'T') out[len - i - 1] = 'A';
    }
    string outString(out);
    free(out);
    return outString;
}

inline double readTimer() {
    return clock() / (double) CLOCKS_PER_SEC;
}

inline string currentDateTime() {
    // Get current date/time, format is YYYY-MM-DD HH:mm:ss
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof (buf), "%Y-%m-%d %X\n", &tstruct);
    return buf;
}


inline int countOutArcs(int node) {
    return (adjList.at(node)).size();
}


inline char boolToCharSign(bool sign) {
    return (sign == true) ? '+' : '-';
}


int maximumUnitigLength(){
    int m = 0;
    for(unitig_struct_t u: unitigs){
        if(u.ln > m){
            m = u.ln;
        }
    }
    return m;
}




int read_unitig_file(const string& unitigFileName, vector<unitig_struct_t>& unitigs) {
    ifstream unitigFile;
    unitigFile.open(unitigFileName);
    
    string line;
    
    int nodeNum;
    char lnline[20];
    char kcline[20];
    char kmline[20];
    char edgesline[100000];
    bool doCont = false;
    
    
    getline(unitigFile, line);
    
    do {
        unitig_struct_t unitig_struct;
        edgesline[0] = '\0';
        sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);
        
        // @@DEBUG
        //        if(unitig_struct.serial == 1241914){
        //            cout<<line<<endl;
        //            cout<<edgesline<<endl;
        //        }
        //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
        
        sscanf(lnline, "%*5c %d", &unitig_struct.ln);
        sscanf(kcline, "%*5c %d", &unitig_struct.kc);
        sscanf(kmline, "%*5c %f", &unitig_struct.km);
        
        char c1, c2;
        stringstream ss(edgesline);
        
        vector<edge_t> edges;
        while (getline(ss, line, ' ')) {
            if (delSpaces(line).length() != 0) {
                if(DBGFLAG==VERIFYINPUT){
                    cout<<line<<endl;
                }
                //                if(unitig_struct.serial == 1241914){
                //                    cout<<line<<endl;
                //                }
                sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                edge_t newEdge;
                
                bool DELSELFLOOP=true;
                if(DELSELFLOOP){
                    if((unitig_struct.serial)!= nodeNum){
                        newEdge.left = charToBool(c1);
                        newEdge.right = charToBool(c2);
                        newEdge.toNode = nodeNum;
                        edges.push_back(newEdge);
                    }
                }else{
                    newEdge.left = charToBool(c1);
                    newEdge.right = charToBool(c2);
                    newEdge.toNode = nodeNum;
                    edges.push_back(newEdge);
                }
                
            }
            
        }
        adjList.push_back(edges);
        
        
        
        doCont = false;
        while (getline(unitigFile, line)) {
            if (line.substr(0, 1).compare(">")) {
                unitig_struct.sequence = unitig_struct.sequence + line;
                unitigs.push_back(unitig_struct);
            } else {
                doCont = true;
                break;
            }
        }
    } while (doCont);
    
    
    unitigFile.close();
    
    cout << "Complete reading input unitig file (bcalm2 file)." << endl;
    return EXIT_SUCCESS;
}


// @@ --- ALL PRINTING CODE --- //
void printBCALMGraph(vector<vector<edge_t> > adjList) {
    for (int i = 0; i < adjList.size(); i++) {
        cout << i << "# ";
        for (edge_t edge : adjList.at(i)) {
            cout << boolToCharSign(edge.left) << ":" << edge.toNode << ":" << boolToCharSign(edge.right) << ", ";
        }
        cout << endl;
    }
}

void printAllBCALMSequences(vector<unitig_struct_t> unitigs) {
    for (unitig_struct_t unitig : unitigs) {
        cout << unitig.serial << ": " << unitig.ln << " " << unitig.sequence.length() << endl;
    }
}


// might be useful for doing some visualization.
bool canReachSinkSource(int v, bool visited[], bool sign)
{
    // Mark the current node as visited and
    // print it
    
    visited[v] = true;
    //cout << v << " ";
    bool reachable = false;
    
    if(global_plusoutdegree[v] == 0 && global_plusindegree[v] != 0){
        //cout<<v<<"is sink.";
        return true;//sink
        
    }
    if(global_plusindegree[v] == 0 && global_plusoutdegree[v] != 0){
        //cout<<v<<"is source.";
        return true;//source
    }
    if(global_indegree[v] == 0){
        //cout<<v<<"is isolated.";
        return true;//isolated
    }
    
    
    // Recur for all the vertices adjacent
    // to this vertex
    vector<edge_t>::iterator i;
    for (i = adjList[v].begin(); i != adjList[v].end(); ++i){
        
        if (!visited[(*i).toNode] && sign==(*i).left){
            reachable = canReachSinkSource((*i).toNode, visited, (*i).right);
            if(reachable==true){
                return true;
            }
        }
        
    }
    return reachable;
    
}


set<int> extractIntegerWords(string str)
{
    set<int> retset;
    stringstream ss;
    
    /* Storing the whole string into string stream */
    ss << str;
    
    /* Running loop till the end of the stream */
    string temp;
    int found;
    while (!ss.eof()) {
        
        /* extracting word by word from stream */
        ss >> temp;
        
        /* Checking the given word is integer or not */
        if (stringstream(temp) >> found)
            retset.insert(found);
        
        /* To save from space at the end of string */
        temp = "";
    }
    return retset;
}

void makeGraphDot(string ipstr){
    FILE * fp;
    
    fp = fopen ("/Users/Sherlock/Downloads/graphviz-2.40.1/graph.gv", "w+");
    
    fprintf(fp, "digraph d {\n");
    //string ipstr = "20 19 18";
    set<int> verticesMarked;
    set<int> vertices = extractIntegerWords(ipstr) ;
    set<int> vMarked(vertices.begin(), vertices.end());
    //set<pair<int, int> > edges;
    
    for(int x: vertices){
        if(x>=adjList.size()){
            cout<<"wrong, do again"<<endl;
            return;
        }
        vector<edge_t> adjX = adjList[x];
        for(edge_t ex: adjX){
            vertices.insert(ex.toNode);
            fprintf(fp, "%d -> %d[taillabel=\"%d\", headlabel=\"%d\", arrowhead=\"none\"]\n", x, ex.toNode, ex.left, !ex.right);
            
            //            pair<int, int> p;
            //            if(x < ex.toNode){
            //                p.first = x;
            //                p.second = ex.toNode;
            //            }else{
            //                p.second = x;
            //                p.first = ex.toNode;
            //            }
            //            edges.insert(p);
        }
    }
    for(int x: vertices){
        if(vMarked.count(x)>0){
            fprintf(fp, "%d [label=\"%d\", color=\"red\"]\n", x, x);
        }else{
            fprintf(fp, "%d [label=\"%d\"]\n", x, x);
        }
        
    }
    
    //for all int in list
    //make a list of neighbors add them
    
    
    fprintf(fp, "}\n");
    
    fclose(fp);
}


ofstream tipFile;
ofstream tipDebugFile;


#endif /* global_h */
