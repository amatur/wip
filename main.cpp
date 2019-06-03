// --- VERSION 1.7 ----
// Bug fixed and validated
// input format input.txt
#include <cmath>
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
using namespace std;

int K = 31;
string UNITIG_FILE = "/Volumes/FAT32/chol31/list_reads.unitigs.fa";
//        K = 31;
//        UNITIG_FILE = "/Users/Sherlock/cse566_2/exclude/staph31/list_reads.unitigs.fa";
//    K = 11;
//    UNITIG_FILE = "/Users/Sherlock/cse566_2/data/list_reads.unitigs.fa";

enum DEBUGFLAG_T { NONE = 0,  UKDEBUG = 0, VERIFYINPUT = 1, INDEGREEPRINT = 2, DFSDEBUGG = 3, PARTICULAR = 4, OLDNEWMAP = 9, PRINTER = 10, SINKSOURCE = 12};

enum ALGOMODE_T { BASIC = 0, INDEGREE_DFS = 1, INDEGREE_DFS_1 = 2, OUTDEGREE_DFS = 3, OUTDEGREE_DFS_1 = 4, INDEGREE_DFS_INVERTED = 5};

DEBUGFLAG_T DBGFLAG = UKDEBUG;
ALGOMODE_T ALGOMODE = OUTDEGREE_DFS;

string mapmode[] = {"random", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs"
};



typedef unsigned char uchar;

typedef struct {
    int serial;
    int startPos;
    int endPos;
} new_node_info_t;

typedef struct {
    int serial;
    string sequence;
    int ln;
    int kc;
    float km;
} unitig_struct_t;

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

int isolated_node_count = 0;
int sink_count = 0;
int source_count = 0;
int onecount = 0;


struct node_sorter {
    int node;
    int sortkey;
    //bool operator() (struct node_sorter  i, struct node_sorter  j) { return (i.sortkey<j.sortkey);}
};
bool sort_by_key (struct node_sorter i, struct node_sorter j) { return (i.sortkey<j.sortkey); }

int* global_indegree;
int* global_outdegree;
int* global_plusindegree;
int* global_plusoutdegree;


vector<vector<edge_t> > adjList;
vector<vector<edge_t> > reverseAdjList;
vector<vector<newEdge_t> > newAdjList;
vector<edge_both_t> resolveLaterEdges;
vector<unitig_struct_t> unitigs;
map<int, string> newSequences;
vector<list<int> > newToOld;


inline string plus_strings(const string& a, const string& b, size_t kmersize) {
    if (a == "") return b;
    if (b == "") return a;
    string ret = a + b.substr(kmersize - 1, b.length() - (kmersize - 1));
    return ret;
}

string delSpaces(string &str) {
    str.erase(std::remove(str.begin(), str.end(), ' '), str.end());
    return str;
}

bool charToBool(char c) {
    if (c == '+') {
        return true;
    } else {
        if (c != '-') cout << "ERRRRRROOR!" << endl;
        return false;
    }
}

string reverseComplement(string base) {
    int len = base.length();
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

double readTimer() {
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

int countInArcs(int node) {
    //system("grep -o -i :2: unitigs.fa | wc -l > incount.txt");
    int count;
    string line;
    string countFile = "incount.txt";
    //string inputFile = "unitigs.fa";
    string inputFile = UNITIG_FILE;
    ostringstream stringStream;
    stringStream << "grep -o -i :" << node << ": " << inputFile << " | wc -l > " << countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    getline(cf, line);
    line = delSpaces(line);
    sscanf(line.c_str(), "%d", &count);
    cf.close();
    return count;
}

int countOutArcs(int node) {
    return (adjList.at(node)).size();
}

int maximumUnitigLength() {
    //grep '>' list_reads.unitigs.fa | cut -d : -f3 | cut -d ' ' -f1 | sort -n | tail -n 1 > incount.txt
    int count;
    string line;
    string countFile = "incount.txt";
    ostringstream stringStream;
    stringStream << "grep '>' " << UNITIG_FILE << " | cut -d : -f3 | cut -d ' ' -f1 | sort -n | tail -n 1 > " << countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    getline(cf, line);
    line = delSpaces(line);
    sscanf(line.c_str(), "%d", &count);
    cf.close();
    return count;
}

inline char boolToCharSign(bool sign) {
    return (sign == true) ? '+' : '-';
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
        cout << unitig.serial << ": " << unitig.ln << " " << unitig.sequence.length() << endl; //sequence only^
        // full print
        //cout<<unitig.serial<<": "<<unitig.ln<<" "<<unitig.sequence.length()<<":"<<unitig.sequence<<endl;
    }
}





class Graph {
public:
    int V = adjList.size();
    int countNewNode;
    int time;
    
    char* color;
    int* p;
    bool* nodeSign;
    new_node_info_t* oldToNew;
    bool* saturated;
    struct node_sorter * indegree;
    struct node_sorter * outdegree;
    
    Graph() {
        color = new char[V];
        p = new int[V];
        nodeSign = new bool[V];
        oldToNew = new new_node_info_t[V];
        saturated = new bool[V];
        indegree = new struct node_sorter[V];
        global_indegree = new int[V];
        global_outdegree = new int[V];
        global_plusindegree = new int[V];
        global_plusoutdegree = new int[V];
        
        for (int i = 0; i < V; i++) {
            oldToNew[i].serial = -1;
            saturated[i] = false;
            indegree[i].sortkey = 0;
            indegree[i].node = i;
            global_indegree[i] = 0;
            global_outdegree[i] = 0;
            global_plusindegree[i] = 0;
            global_plusoutdegree[i] = 0;
        }
    }
    
    void indegreePopulate(){
        int x = 0;
        for(vector<edge_t> elist: adjList){
            for(edge_t e: elist){
                global_indegree[e.toNode] += 1;
                indegree[e.toNode].sortkey = indegree[e.toNode].sortkey + 1;
                if(e.right == true){
                    global_plusindegree[e.toNode] += 1;
                }
                if(e.left == true){
                    global_plusoutdegree[x] += 1;
                }
                
            }
            global_outdegree[x] = elist.size();
            x++;
        }
        
        for(int i = 0; i<V; i++){
            int minusindegree = (global_indegree[i] - global_plusindegree[i] );
            int minusoutdegree = (global_outdegree[i] - global_plusoutdegree[i] );
            
            if(DBGFLAG == SINKSOURCE){
                cout<<i<<"is ";
            }
//            if(global_plusindegree[i] != 0 && global_plusoutdegree[i] == 0){
//                sink_count++;
//                cout<<"sink, ";
//            }else if(minusindegree != 0 && minusoutdegree == 0){
//                sink_count++;
//                cout<<"sink, ";
//            }
//
//            if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
//                source_count++;
//                cout<<"source, ";
//            }else if(minusindegree == 0 && minusoutdegree != 0){
//                source_count++;
//                cout<<"source, ";
//            }
            
            
            if(global_plusoutdegree[i] == 0){
                sink_count++;
                if(DBGFLAG == SINKSOURCE){
                    cout<<"sink, ";
                }
                
            }else if(minusoutdegree == 0){
                sink_count++;
                if(DBGFLAG == SINKSOURCE){
                    cout<<"sink, ";
                }
            }
            
            if(global_plusindegree[i] == 0){
                source_count++;
                if(DBGFLAG == SINKSOURCE){
                    cout<<"source, ";
                }
            }else if(minusindegree == 0){
                source_count++;
                if(DBGFLAG == SINKSOURCE){
                    cout<<"source, ";
                }
            }
            
            
            global_outdegree[i] += global_indegree[i];
            if(global_indegree[i] == 0){
                isolated_node_count++;
                if(DBGFLAG == SINKSOURCE){
                    cout<<"isolated, ";
                }
            }
            if(global_indegree[i] == 1){
                onecount++;
            }
            
            if(DBGFLAG == SINKSOURCE){
                cout<<endl;
            }
            
        }
    }
    
    
    void DFS_visit(int u) {
        stack<edge_t> s;
        edge_t uEdge;
        uEdge.toNode = u;
        s.push(uEdge);
        
        while (!s.empty()) {
            edge_t xEdge = s.top();
            
            
            
            int x = xEdge.toNode;
            s.pop();
            
            
            
            
            if (color[x] == 'w') {
                //Original DFS code
                time = time + 1;
                color[x] = 'g';
                s.push(xEdge);
                vector<edge_t> adjx = adjList.at(x);
                
                
                if(ALGOMODE == INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] < global_indegree[rhs.toNode];
                         });
                }
                if(ALGOMODE == INDEGREE_DFS_INVERTED){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] > global_indegree[rhs.toNode];
                         });
                }
                if (ALGOMODE == OUTDEGREE_DFS){
                    if(p[x] == -1){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_outdegree[lhs.toNode] < global_outdegree[rhs.toNode];
                             });
                    }else if (nodeSign[x] == false){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_outdegree[lhs.toNode] - global_plusoutdegree[lhs.toNode] < global_outdegree[rhs.toNode] - global_plusoutdegree[rhs.toNode];
                             });
                    }else if (nodeSign[x] == true){
                        sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                             {
                                 return global_plusoutdegree[lhs.toNode] < global_plusoutdegree[rhs.toNode];
                             });
                    }
                }
                
                
                // Now our branching code ::
                
                // For a white x
                // Consider 2 case:
                // Case 1. p[x] = -1, it can happen in two way, x is the first one ever in this connected component, or no one wanted to take x
                // either way, if p[x] = -1, i can be representative of a new node in new graph
                // Case 2. p[x] != -1, so x won't be the representative/head of a newHome. x just gets added to its parent's newHome.
                int u = unitigs.at(x).ln; //unitig length
                
                if (p[x] == -1) {
                    
                    list<int> xxx;
                    xxx.push_back(x);
                    newToOld.push_back(xxx);
                    oldToNew[x].serial = countNewNode++; // countNewNode starts at 0, then keeps increasing
                    
                    //make the sequence
                    //NOT CORRECT? I am not sure
                    if(nodeSign[x]==false){
                        newSequences[oldToNew[x].serial] = reverseComplement(unitigs.at(x).sequence);
                    }else{
                        newSequences[oldToNew[x].serial] = (unitigs.at(x).sequence);
                    }
                    
                    
                    oldToNew[x].startPos = 1;
                    if (u < K) {
                        oldToNew[x].endPos = 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPos = u - K + 1;
                    }
                    
                } else {
                    
                    newToOld[oldToNew[p[x]].serial].push_back(x);
                    oldToNew[x].serial = oldToNew[p[x]].serial;
                    oldToNew[x].startPos = oldToNew[p[x]].endPos + 1;
                    if (u < K) {
                        oldToNew[x].endPos = oldToNew[x].startPos + 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPos = u - K + (oldToNew[x].startPos); //check correctness
                    }
                    
                    // x says: Now that I know where my newHome is: I can extend my parent's sequence
                    // Is it more complicated than this?
                    string parentSeq = newSequences[oldToNew[x].serial];
                    string childSeq = unitigs.at(x).sequence;
                    
                    // Is it CORRECT? just for testing now
                    if(nodeSign[x]==false){
                        childSeq = reverseComplement(childSeq);
                    }
                    newSequences[oldToNew[x].serial] = plus_strings(parentSeq, childSeq, K);
                }
                
                // x->y is the edge, x is the parent we are extending
                for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                    int y = yEdge.toNode;
                    
                    if (DBGFLAG == DFSDEBUGG) {
                        cout << "Edge " << x << "->" << y << endl;
                    }
                    
                    //Normal DFS
                    if (color[y] == 'w') {
                        s.push(yEdge);
                    }
                    
                    if(DBGFLAG == PARTICULAR){
                        // DEBUGGGING a particular edge
                        if (y == 2 && x == 0) {
                            cout << "Saturated? " << saturated[x] << endl;
                        }
                    }
                    
                    
                    
                    //handle self-loop, self-loop will always be an extra edge
                    if (y == x) {
                        edge_both_t e;
                        e.edge = yEdge;
                        e.fromNode = x;
                        resolveLaterEdges.push_back(e);
                    } else if (saturated[x]) {
                        // Since x is saturated, we only add resolveLater edges
                        // no need to check for consistency
                        if (y != p[x]) {
                            edge_both_t e;
                            e.edge = yEdge;
                            e.fromNode = x;
                            resolveLaterEdges.push_back(e);
                        }
                    } else {
                        // If x has space to take a child, meaning x is not saturated
                        // hunting for potential child
                        
                        if (color[y] == 'w' && p[y] == -1) {
                            // y has white color & got no parent => means it's homeless, so let's see if we can take it as a child of x
                            //But just see if it is eligible to be a child, i.e. is it consistent (sign check)?
                            
                            //2 case, Does x's child have grandparent?
                            // No.
                            if (p[x] == -1) {
                                // case 1: child has no grandparent
                                // so extend path without checking any sign
                                nodeSign[x] = yEdge.left;
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true; //found a child
                                
                                
                                //TESTED NOT YET
                                //                                if (nodeSign[y] == false) {
                                //                                    unitigs.at(y).sequence = reverseComplement(unitigs.at(y).sequence);
                                //                                }
                                
                                //Yes.
                            } else if (nodeSign[x] == yEdge.left) {
                                // case 2: child (=y) has grandparent, i.e. x's parent exists
                                nodeSign[y] = yEdge.right;
                                p[y] = x;
                                saturated[x] = true; //found a child
                                
                                //TESTED NOT YET
                                //                                if (nodeSign[y] == false) {
                                //                                    unitigs.at(y).sequence = reverseComplement(unitigs.at(y).sequence);
                                //                                }
                                
                            } else {
                                // do we reach this case?
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                                resolveLaterEdges.push_back(e);
                            }
                            
                        } else {
                            if (y != p[x]) {
                                edge_both_t e;
                                e.edge = yEdge;
                                e.fromNode = x;
                                resolveLaterEdges.push_back(e);
                                if (DBGFLAG == PARTICULAR) {
                                    // DEBUGGGING a particular edge
                                    if (y == 2 && x == 0) {
                                        cout << "Saturated? " << saturated[x] << endl;
                                    }
                                }
                            } else {
                                if (DBGFLAG == PARTICULAR) {
                                    // DEBUGGGING a particular edge
                                    if (y == 2 && x == 0) {
                                        cout << "Saturated? " << saturated[x] << endl;
                                    }
                                }
                            }
                            
                            
                        }
                    }
                }
                
            } else if (color[x] == 'g') {
                time = time + 1;
                color[x] = 'b';
            }
            
        }
    }
   
    
    void DFS() {
        indegreePopulate();
        if (ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1 ){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            copy(myvector.begin(), myvector.end(), indegree);
            
            if(DBGFLAG == INDEGREEPRINT){
                cout<<"print in degrees"<<endl;
                for(int i = 0; i<V; i++){
                    cout<<indegree[i].node<<"->"<<indegree[i].sortkey<<endl;
                }
            }
        }
        
        if (ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_outdegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            copy(myvector.begin(), myvector.end(), indegree);
        }
        
        
        for (int i = 0; i < V; i++) {
            color[i] = 'w';
            p[i] = -1;
        }
        
        for (int j = 0; j < V; j++) {
            int i;
            if(ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1 || ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1){
                i = indegree[j].node;
            }else{
                i = j;
            }
            
            if (color[i] == 'w') {
                if(DBGFLAG == DFSDEBUGG){
                    cout<<"visit start of node: "<<i<<endl;
                }
                DFS_visit(i);
            }
        }
        
        
        //reuse the colors
        for (int i = 0; i < V; i++) {
            color[i] = 'w';
        }
        
        
        for (int i = 0; i < countNewNode; i++) {
            newAdjList.push_back(vector<newEdge_t>());
        }
        
        for (edge_both_t e : resolveLaterEdges) {
            // e.start -> e.end : they have different newHome
            // look at the sign of both end points. if from sign and edge from are not equal, then revert the edge label
            // We consider that ACT and CTT has an edge.
            // ACT(+)  (+)->(-)  AAG(+) : this is fine, just add the edge
            // ACT(-)  (+)->(-)  AAG(+) : you need to convert it, because => AGT  (+)->(-)  AAG => this is not correct:
            // to fix this AGT (-)->(-) AAG
            // don't touch the node sign: just fix the edge signs
            
            int x = e.fromNode;
            int u = unitigs.at(x).ln;
            newEdge_t newEdge;
            newEdge.kmerEndIndex = oldToNew[e.edge.toNode].startPos;
            newEdge.kmerStartIndex = oldToNew[x].endPos;
            
            if(e.fromNode!=e.edge.toNode) {
                if (nodeSign[e.fromNode] != e.edge.left) {
                    e.edge.left = !e.edge.left;
                    newEdge.kmerStartIndex = oldToNew[x].startPos;
                }
                if (nodeSign[e.edge.toNode] != e.edge.right) {
                    e.edge.right = !e.edge.right;
                    newEdge.kmerEndIndex = oldToNew[e.edge.toNode].endPos;
                }
            }
            
            
            newEdge.edge = e.edge;
            newEdge.edge.toNode = oldToNew[newEdge.edge.toNode].serial;
            
            newAdjList[oldToNew[x].serial].push_back(newEdge);
            
            if(DBGFLAG == OLDNEWMAP){
                cout << "old: " << x << "->" << e.edge.toNode << ", new:" << " (" << oldToNew[x].serial << "->" << newEdge.edge.toNode << ")" << endl;
                
            }
        }
    }
    
    
    ~Graph() {
        delete [] color;
        delete [] p;
        delete [] nodeSign;
        delete [] oldToNew;
        delete [] saturated;
        delete [] indegree;
        delete [] global_indegree;
        delete [] global_outdegree;
        delete [] global_plusindegree;
        delete [] global_plusoutdegree;
    }
};


void printNewGraph(Graph &G){
    list<int> *newToOld = new list<int>[G.countNewNode];
    
    // PRINT NEW GRAPH
    for (int i = 0; i < G.V; i++) {
        newToOld[G.oldToNew[i].serial].push_back(i);
        //cout << "old " << i << "-> new" << G.oldToNew[i].serial << endl;
    }
    
    for (int i = 0; i < G.countNewNode; i++) {
        list<int> adj = newToOld[i];
        
        cout<<"new " << i<<": old (index) ";
        for(int val : adj){
            cout<<val<<" ";
        }
        cout<<endl;
        
    }
    
    delete [] newToOld;
    
}

void formattedOutput(Graph &G){
    
    string stitchedUnitigs = "stitchedUnitigs.txt";
    ofstream myfile;
    myfile.open (stitchedUnitigs);
    //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
        myfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
        vector<newEdge_t> edges = newAdjList.at(newNodeNum);
        for(newEdge_t edge: edges){
            myfile<<"L:" << edge.kmerStartIndex << ":" << boolToCharSign(edge.edge.left) << ":" << edge.edge.toNode << ":" << boolToCharSign(edge.edge.left) <<":"<< edge.kmerEndIndex <<" ";
        }
        myfile<<endl;
        myfile<<newSequences[newNodeNum];
        myfile<<endl;
    }
    //myfile << '>' << newNodeNum <<">0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- " ;
    myfile.close();
    
}

int getNodeNumFromFile(string filename) {
    //grep '>' list_reads.unitigs.fa | tail -n 1
    int nodeCount;
    string countFile = "incount.txt";
    ostringstream stringStream;
    stringStream << "grep '>' " << filename << " | tail -n 1" << countFile;
    string copyOfStr = stringStream.str();
    system(copyOfStr.c_str());
    ifstream cf;
    cf.open(countFile);
    string line;
    getline(cf, line);
    sscanf(line.c_str(), "%*c %d", &nodeCount);
    cf.close();
    return nodeCount;
}

int get_data(const string& unitigFileName,
             uchar*& data,
             vector<unitig_struct_t>& unitigs,
             uint64_t& char_count
             ) {
    ifstream unitigFile;
    try {
        unitigFile.open(unitigFileName);
        //            if(unitigFile==NULL){
        //                throw "ERROR: File does not exist!!!";
        //            }else{
        //                cout<<"File opened successfully."<<endl;
        //            }
        
    } catch (const char* msg) {
        cout << msg << endl;
    }
    
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
        
        sscanf(line.c_str(), "%*c %d %s  %s  %s %[^\n]s", &unitig_struct.serial, lnline, kcline, kmline, edgesline);
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
                sscanf(line.c_str(), "%*2c %c %*c %d  %*c  %c", &c1, &nodeNum, &c2); //L:-:0:-
                edge_t newEdge;
                
                newEdge.left = charToBool(c1);
                newEdge.right = charToBool(c2);
                newEdge.toNode = nodeNum;
                edges.push_back(newEdge);
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
    
    cout << "Complete reading input." << endl;
    return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
    
//    string line;
//    ifstream afile ("input.txt");
//    if (afile.is_open())
//    {
//        getline (afile, UNITIG_FILE);
//        getline (afile, line);
//        K = stoi(line);
//        afile.close();
//    }
    
    
    uint64_t char_count;
    uchar *data = NULL;
    
    
    double startTime = readTimer();
    
    cout << "Starting reading file: " << UNITIG_FILE << endl;
    if (EXIT_FAILURE == get_data(UNITIG_FILE, data, unitigs, char_count)) {
        return EXIT_FAILURE;
    }
    cout<<K<<endl;
    //printBCALMGraph(adjList);
    
    double TIME_READ_SEC = readTimer() - startTime;
    
    Graph G;
    G.DFS();
    
    
    if(DBGFLAG == PRINTER){
        printBCALMGraph(adjList);
        printNewGraph(G);
        
        for(int i = 0; i< G.countNewNode; i++){
            cout<<"new ->" <<i<<" ";
            for(int x: newToOld[i]){
                cout<<x<<" ";
            }
            cout<<endl;
        }
    }
    
    
    //fix sequences
    for(int i = 0; i< G.countNewNode; i++){
        string s = "";
        for(int x: newToOld[i]){
            if(G.nodeSign[x] == false){
                s = plus_strings(s, reverseComplement(unitigs.at(x).sequence), K);
            }else{
                s = plus_strings(s, (unitigs.at(x).sequence), K);
            }
        }
        newSequences[i] = s;
        //cout<<endl;
    }
    
    
    // COLLECT STATISTICS
    
    //count total number of edges
    int E = 0;
    for (int i = 0; i < G.V; i++) {
        E += adjList[i].size();
    }
    
    int E_new = resolveLaterEdges.size();
    int V = G.V;
    int V_new = G.countNewNode;
    
    int C = 0;
    int C_new = 0;
    for (unitig_struct_t unitig : unitigs) {
        C += unitig.ln;
    }

    
    
    map<int, string>::iterator it;
    int maxlen = 0;
    for (it = newSequences.begin(); it != newSequences.end(); it++)
    {
        C_new += (it->second).length();
        if((it->second).length() >maxlen){
            maxlen =(it->second).length();
        }
        
    }
    

    
    double TIME_TOTAL_SEC = readTimer() - startTime;
    
    
    // For collecting stats
    int U_MAX = maximumUnitigLength();
    int EDGE_INT_DTYPE_SIZE;
    if (U_MAX - K + 1 > 0) {
        EDGE_INT_DTYPE_SIZE = log2(U_MAX - K + 1);
    } else {
        EDGE_INT_DTYPE_SIZE = log2(K);
    }
    EDGE_INT_DTYPE_SIZE = ceil(EDGE_INT_DTYPE_SIZE / 8.0);
    
    
    int ACGT_DTYPE_SIZE = 1; // 1 byte to store each char
    int NODENUM_DTYPE_SIZE = 8; // 1 byte to store each char
    int SIGN_DTYPE_SIZE = 1; // 1 byte for sign information (+,- can be stored in 1 byte)
    int spaceBefore = C * ACGT_DTYPE_SIZE + E * (NODENUM_DTYPE_SIZE + SIGN_DTYPE_SIZE);
    int save = (C - C_new) * ACGT_DTYPE_SIZE + (E - E_new)*(NODENUM_DTYPE_SIZE + SIGN_DTYPE_SIZE);
    int overhead = (E_new)*(2 * EDGE_INT_DTYPE_SIZE);
    float persaved = ((save - overhead)*1.0 / spaceBefore) * 100.0;
    float upperbound = (1-((C-(K-1)*(G.V - max(sink_count, source_count)*1.0))/C))*100.0;
    float saved_c = (1-(C_new*1.0/C))*100.0;
    
    formattedOutput(G);
    printf("%d \t %d \t %d \t %d \t %d \t %d \t %f \t %f \t %.2f%% \t %d \t %d \t %d \t %f \t %f \t %d \t %d \t %d \t %d \t %.2f%% \t %.2f%% \t %s\n", V, V_new, E, E_new, C, C_new, spaceBefore / 1024.0, (save - overhead) / 1024.0, persaved, U_MAX, maxlen, K, TIME_READ_SEC, TIME_TOTAL_SEC, isolated_node_count, onecount, sink_count, source_count, upperbound, saved_c, mapmode[ALGOMODE].c_str());
    
    return EXIT_SUCCESS;
}
