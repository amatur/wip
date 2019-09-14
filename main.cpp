// --- VERSION 5.2 ----
// - Aug 28
// forward ext + two way + bracket error

//Caution:
//removed all self-loops


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


using namespace std;

typedef tuple<int,int,int, int> mytuple; // uid, walkid, pos, isTip

bool sort_by_walkId (const mytuple &lhs, const mytuple &rhs){
    return get<1>(lhs) < get<1>(rhs);
}
bool sort_by_pos (const mytuple &lhs, const mytuple &rhs){
    return get<2>(lhs) < get<2>(rhs);
}
bool sort_by_tipstatus (const mytuple &lhs, const mytuple &rhs){
    return get<3>(lhs) < get<3>(rhs);
}

int K = 55;
string UNITIG_FILE = "/Volumes/exFAT/data2019/chol/55/list_reads.unitigs.fa";
int C_twoway = 0;
int C_bracketed = 0;
int V_twoway = 0;
int V_new = 0;
int C_new = 0;
int V_bracketed = 0;


enum DEBUGFLAG_T { NONE = 0,  UKDEBUG = 0, VERIFYINPUT = 1, INDEGREEPRINT = 2, DFSDEBUGG = 3, PARTICULAR = 4, NODENUMBER_DBG = 5, OLDNEWMAP = 9, PRINTER = 10, SINKSOURCE = 12};

enum ALGOMODE_T { BASIC = 0, INDEGREE_DFS = 1, INDEGREE_DFS_1 = 2, OUTDEGREE_DFS = 3, OUTDEGREE_DFS_1 = 4, INDEGREE_DFS_INVERTED = 5, PLUS_INDEGREE_DFS = 6, RANDOM_DFS = 7, NODEASSIGN = 8, SOURCEFIRST = 9, TWOWAYEXT = 10, PROFILE_ONLY = 11, EPPRIOR=12, GRAPHPRINT = 13, TIGHTUB = 14, BRACKETCOMP = 15};

bool FLG_NEWUB = true;

DEBUGFLAG_T DBGFLAG = NONE; //NODENUMBER_DBG
ALGOMODE_T ALGOMODE = TWOWAYEXT;

string mapmode[] = {"basic", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "two_way_extension", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "bracket_comp"
};
string modefilename[] = {"Fwd", "indegree_dfs", "indegree_dfs_initial_sort_only", "outdegree_dfs", "outdegree_dfs_initial_sort_only", "inverted_indegree_dfs", "plus_indegree_dfs", "random_dfs", "node_assign", "source_first", "", "profile_only", "endpoint_priority", "graph_print", "tight_ub", "Tip"
};

typedef unsigned char uchar;



typedef struct {
    int serial = -1;
    int startPosWithKOverlap;
    int endPosWithKOVerlap;
    bool isWalkEnd = false;
    int pos_in_walk = -1;
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
int sharedparent_count = 0;
int sharparentCntRefined = 0;
int onecount = 0;


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

map<pair <int, int>, int> inOutCombo;

vector<vector<edge_t> > adjList;
vector<vector<edge_t> > reverseAdjList;

vector<edge_both_t> resolveLaterEdges;
//vector<edge_both_t> sinkSrcEdges;
vector<unitig_struct_t> unitigs;
map<int, string> newSequences;
map<int, string> newNewSequences; //int is the unitig id (old id)
set<int> newNewMarker;

vector<list<int> > newToOld;
vector<int> walkFirstNode; //given a walk id, what's the first node of that walk
unordered_map<int, vector<edge_t> > sinkSrcEdges; //int is the unitig id (old id)

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


int countOutArcs(int node) {
    return (adjList.at(node)).size();
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
        cout << unitig.serial << ": " << unitig.ln << " " << unitig.sequence.length() << endl;
    }
}



class GroupMerger {
public:
    map<int, bool> fwdVisited;
    map<int, bool> bwdVisited;
    map<int, int> bwdWalkId;
    map<int, int> fwdWalkId;
    GroupMerger() {
    }
    void connectGroups(int from, int to){
        fwdVisited[from] = false;
        bwdVisited[to] = false;
        fwdWalkId[from] = to;
        bwdWalkId[to] = from;
    }
    ~GroupMerger() {
    }
};


class DisjointSet {
    unordered_map<int, int> parent;
    
public:
    DisjointSet() {
    }
    void make_set(int id) {
        this->parent[id] = -1;
    }
    
    void Union(int xId, int yId) {
        int xset = find_set(xId);
        int yset = find_set(yId);
        
        if(xset != yset)
        {
            parent[xset] = yset;
        }
    }
    
    int find_set(int id) {
        if (parent[id] == -1)
            return id;
        return find_set(parent[id]);
    }
    ~DisjointSet(){
    }
    
};


class Graph {
public:
    size_t V = adjList.size();
    int countNewNode = 0;
    int time = 0;
    
    char* color;
    int* p;
    bool* nodeSign;
    new_node_info_t* oldToNew;
    bool* saturated;
    struct node_sorter * indegree;
    struct node_sorter * outdegree;
    bool* countedForLowerBound;
    DisjointSet disSet;
    GroupMerger gmerge;
    
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
        global_issinksource = new int[V];
        global_priority = new int[V];
        countedForLowerBound = new bool[V];
        
        for (int i = 0; i < V; i++) {
            
            if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP ){
                disSet.make_set(i);
            }
            
            oldToNew[i].serial = -1;
            saturated[i] = false;
            indegree[i].sortkey = 0;
            indegree[i].node = i;
            global_indegree[i] = 0;
            global_outdegree[i] = 0;
            global_plusindegree[i] = 0;
            global_plusoutdegree[i] = 0;
            global_issinksource[i] = 0;
            global_priority[i] = 0;
            countedForLowerBound[i] = false;
        }
    }
    
    inline bool sRight(edge_t plusminusedge){
        return !(plusminusedge.right == true);
    }
    
    inline bool sLeft(edge_t plusminusedge){
        return (plusminusedge.left == true);
    }
    
    void indegreePopulate(){
        int xc = 0;
        for(vector<edge_t> elist: adjList){
            for(edge_t e: elist){
                global_indegree[e.toNode] += 1;
                indegree[e.toNode].sortkey = indegree[e.toNode].sortkey + 1;
                if(e.right == true){
                    global_plusindegree[e.toNode] += 1;
                }
                if(e.left == true){
                    global_plusoutdegree[xc] += 1;
                }
                
            }
            global_outdegree[xc] = elist.size();
            xc++;
        }
        
        for(int i = 0; i<5; i++){
            for(int j = 0; j<5; j++){
                inOutCombo[make_pair(i,j)] = 0;
            }
        }
        for(int i = 0; i<V; i++){
            pair<int, int> a;
            a = make_pair(global_plusindegree[i], global_plusoutdegree[i]);
            inOutCombo[a] = (inOutCombo.count(a)  ? inOutCombo[a] + 1 : 1  );
            
            
            if(DBGFLAG == SINKSOURCE){
                cout<<i<<"is ";
            }
            
            if(global_plusoutdegree[i] == 0 && global_plusindegree[i] != 0){
                sink_count++;
                global_issinksource[i] = 1;
                global_priority[i] = 5;
                countedForLowerBound[i] = true;
                
                if(DBGFLAG == SINKSOURCE){
                    cout<<"sink, ";
                }
                
            }
            if(global_plusindegree[i] == 0 && global_plusoutdegree[i] != 0){
                source_count++;
                global_issinksource[i] = 1;
                global_priority[i] = 5;
                countedForLowerBound[i] = true;
                
                if(DBGFLAG == SINKSOURCE){
                    cout<<"source, ";
                }
            }
            
            if(global_indegree[i] == 0){
                global_issinksource[i] = 1;
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
        
        xc = 0; // current vertex while traversing the adjacency list
        for(vector<edge_t> elist: adjList){
            int neighborCount = 0;
            stack<int> countedNodes;
            set<pair<int, bool> > countedSides;
            //if(true){
            if(FLG_NEWUB){
                //ENDPOINT SIDE UPPER BOUND
     
                for(edge_t e_xy: elist){
                    int y = e_xy.toNode;
                    vector<edge_t> adjY = adjList[y];
                    bool eligible = true;
                    pair<int, bool> pairr;
                    for(edge_t e_from_y : adjY){
                        pairr =make_pair(e_from_y.toNode, sRight(e_xy) );
                        if(e_from_y.toNode!=xc){
                            
                            if(sRight(e_xy) == sLeft(e_from_y)){
                                eligible = false;
                                break;
                            }
                            
                        }
                        
                    }
                    
                    if(eligible){
                        neighborCount++;
                        //if(countedSides.count(pairr)==0)
                        
                        
                        //countedSides.insert(pairr);
                    }
                }
                
                if(global_issinksource[xc] == 1){
                    if(neighborCount>1){
                        sharedparent_count += neighborCount - 1 ;
                    }
                }else{
                    if(neighborCount>2){
                        sharedparent_count += neighborCount - 2 ;
                    }
                }
            }
            //sharedparent_count_wrong =sharedparent_count;
            
            //if(true){
            if(!FLG_NEWUB){
                // OLDER UPPER BOUND CALC
                
                int neighborCount = 0;
                for(edge_t e_xy: elist){
                    int y = e_xy.toNode;
                    
                    //                if(!counted[y] && global_plusindegree[y] == 1 && global_plusoutdegree[y] == 1 && e_xy.right == true ){
                    //                    neighborCnt++;
                    //                    counted[y] = true;
                    //                }
                    
                    if(!countedForLowerBound[y]){
                        vector<edge_t> adjY = adjList[y];
                        bool eligible = true;
                        for(edge_t e_from_y : adjY){
                            if(e_from_y.toNode!=xc){
                                if(sRight(e_xy) == sLeft(e_from_y) ){
                                    eligible = false;
                                    break;
                                }
                            }
                        }
                        if(eligible){
                            countedForLowerBound[y] = true;
                            global_priority[y] = 4;
                            neighborCount++;
                            countedNodes.push(y);
                        }
                        
                    }
                }
                
                
                if(global_issinksource[xc] == 1){
                    if(neighborCount>1){
                        sharedparent_count += neighborCount - 1 ;
                    }else{
                        while(!countedNodes.empty()){
                            countedForLowerBound[countedNodes.top()] = false;
                            countedNodes.pop();
                        }
                    }
                }else{
                    if(neighborCount>2){
                        sharedparent_count += neighborCount - 2 ;
                    }else{
                        while(!countedNodes.empty()){
                            countedForLowerBound[countedNodes.top()] = false;
                            countedNodes.pop();
                        }
                    }
                }
            }
            
            xc++;
        }
    }
    
    
    void DFS_visit(int u) {
        if(ALGOMODE == BRACKETCOMP){
            
            if(global_issinksource[u]==1){
                vector<edge_t> adju = adjList.at(u);
                vector<edge_t> myvector;
                for (edge_t e : adju) {
                    myvector.push_back(e);
                }
                sinkSrcEdges[u] = myvector;
                return;
            }
        }
        
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
                if(ALGOMODE == RANDOM_DFS){
                    random_shuffle ( adjx.begin(), adjx.end() );
                }
                
                if(ALGOMODE == EPPRIOR){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_priority[lhs.toNode]   <  global_priority[rhs.toNode]  ;
                         });
                }
                
                if(ALGOMODE == INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode] < global_indegree[rhs.toNode];
                         });
                }
                
                if(ALGOMODE == PLUS_INDEGREE_DFS){
                    sort( adjx.begin( ), adjx.end( ), [ ]( const edge_t& lhs, const edge_t& rhs )
                         {
                             return global_indegree[lhs.toNode]  - global_plusindegree[lhs.toNode] > global_indegree[lhs.toNode]  - global_plusindegree[rhs.toNode];
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
                    oldToNew[x].finalWalkId = oldToNew[x].serial;
                    
                    
                    //added while doing bracket comp
                    walkFirstNode.push_back(x);
                    
                    
                    //make the sequence
                    //NOT CORRECT? I am not sure
                    if(nodeSign[x]==false){
                        newSequences[oldToNew[x].serial] = reverseComplement(unitigs.at(x).sequence);
                        
                        
                        //newNewSequences[x] = reverseComplement(unitigs.at(x).sequence);
                    }else{
                        newSequences[oldToNew[x].serial] = (unitigs.at(x).sequence);
                        
                        //newNewSequences[x] = (unitigs.at(x).sequence);
                    }
                    
                    oldToNew[x].pos_in_walk = 1;
                    oldToNew[x].startPosWithKOverlap = 1;
                    if (u < K) {
                        oldToNew[x].endPosWithKOVerlap = 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPosWithKOVerlap = u - K + 1;
                    }
                    
                } else {
                    
                    newToOld[oldToNew[p[x]].serial].push_back(x);
                    oldToNew[x].serial = oldToNew[p[x]].serial;
                    oldToNew[x].finalWalkId = oldToNew[x].serial;
                    
                    
                    if(ALGOMODE==TWOWAYEXT || ALGOMODE==BRACKETCOMP ){
                        disSet.Union(x, p[x]);
                    }
                    
                    
                    oldToNew[x].startPosWithKOverlap = oldToNew[p[x]].endPosWithKOVerlap + 1;
                    oldToNew[x].pos_in_walk = oldToNew[p[x]].pos_in_walk + 1;
                    if (u < K) {
                        oldToNew[x].endPosWithKOVerlap = oldToNew[x].startPosWithKOverlap + 1; // do we actually see this? yes
                        if(DBGFLAG == UKDEBUG){
                            cout<< "node: "<< x<<"u< k ***** u = "<<u<<endl;
                        }
                    } else {
                        oldToNew[x].endPosWithKOVerlap = u - K + (oldToNew[x].startPosWithKOverlap); //check correctness
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
                    newNewSequences[x] = newSequences[oldToNew[x].serial] ;
                }
                
                // x->y is the edge, x is the parent we are extending
                for (edge_t yEdge : adjx) { //edge_t yEdge = adjx.at(i);
                    int y = yEdge.toNode;
                    
                    if(ALGOMODE == BRACKETCOMP){
                        if(global_issinksource[y] == true){
                            continue;
                        }
                        
                    }
                    
                    // cout << "Edge " << x << "->" << y << " "<<global_indegree[y] << endl;
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
                            if (p[x] == -1 && ALGOMODE != NODEASSIGN) {
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
                            
                            //merger
                            if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP){
                                // y is not white
                                bool consistentEdge = (nodeSign[y] == yEdge.right && (p[x]==-1 || (p[x]!=-1&& nodeSign[x] == yEdge.left)) );
                                if(p[y]==-1 && consistentEdge && oldToNew[x].serial != oldToNew[y].serial){
                                    
                                    //cout<<"x: "<<x<<":" <<disSet.find_set(x)<<" ";
                                    //cout<<"y: "<<y<<":" <<disSet.find_set(y) <<endl;
                                    
                                    //not in same group already, prevent cycle
                                    if(disSet.find_set(x)!=disSet.find_set(y)){
                                        
                                        nodeSign[x] = yEdge.left;
                                        nodeSign[y] = yEdge.right;
                                        p[y] = x;
                                        saturated[x] = true; //found a child
                                        // oldToNew[y].serial
                                        
                                        disSet.Union(x, y);
                                        //cout<<"x: "<<disSet.find_set(x);
                                        //cout<<"y: "<<disSet.find_set(y);
                                        //cout<<endl;
                                        
                                        //cout<<oldToNew[x].serial<<"+"<<oldToNew[y].serial<<endl;
                                        gmerge.connectGroups(oldToNew[x].serial,oldToNew[y].serial );
                                        
                                    }
                                    
                                }
                                
                            }
                            
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
    
    void DFS() {

        if(ALGOMODE == NODEASSIGN){
            for (int i=0; i<V; i++) {
                nodeSign[i] = true;
                if(global_plusindegree[i]< global_indegree[i] - global_plusindegree[i]){
                    nodeSign[i] = true;
                }
            }
        }
        
        if(ALGOMODE == SOURCEFIRST){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_issinksource[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            copy(myvector.begin(), myvector.end(), indegree);
        }
        
        if(ALGOMODE == EPPRIOR){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_priority[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            copy(myvector.begin(), myvector.end(), indegree);
        }
        
        if(ALGOMODE == INDEGREE_DFS_INVERTED){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key_inverted);
            //random_shuffle ( myvector.begin(), myvector.end() );
            copy(myvector.begin(), myvector.end(), indegree);
            
        }
        
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
        
        
        if(ALGOMODE == RANDOM_DFS){
            for (int i = 0; i < V; i++) {
                indegree[i].node = i;
                indegree[i].sortkey = global_indegree[i];
            }
            vector<struct node_sorter> myvector (indegree, indegree+V);
            sort (myvector.begin(), myvector.end(), sort_by_key);
            random_shuffle ( myvector.begin(), myvector.end() );
            copy(myvector.begin(), myvector.end(), indegree);
            
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
        
        double time_a = readTimer();
        for (int i = 0; i < V; i++) {
            color[i] = 'w';
            p[i] = -1;
        }
        cout<<"Basic V loop time: "<<readTimer() - time_a<<" sec"<<endl;
        
        
        time_a = readTimer();
        for (int j = 0; j < V; j++) {
            int i;
            if(ALGOMODE == OUTDEGREE_DFS || ALGOMODE == OUTDEGREE_DFS_1 || ALGOMODE == INDEGREE_DFS || ALGOMODE == INDEGREE_DFS_1 || ALGOMODE == SOURCEFIRST){
                i = indegree[j].node;
            }else{
                i = j;
            }
            
            
            
            if (color[i] == 'w') {
                if(DBGFLAG == DFSDEBUGG ){
                    cout<<"visit start of node: "<<i<<endl;
                }
                DFS_visit(i);
            }
        }
        cout<<"DFS time: "<<readTimer() - time_a<<" sec"<<endl;
        

        
        
        cout<<"## START stitching strings: "<<endl;
        time_a = readTimer();
        
        //fix sequences
        for(int i = 0; i< countNewNode; i++){
            string s = "";
            for(int x: newToOld[i]){
                if(color[x] != 'l' && color[x] != 'r' ){ // useless in normal code
                    if(nodeSign[x] == false){
                        s = plus_strings(s, reverseComplement(unitigs.at(x).sequence), K);
                    }else{
                        s = plus_strings(s, (unitigs.at(x).sequence), K);
                    }
                }
            }
            
            newSequences[i] = s;
            
            for(int x: newToOld[i]){
                newNewSequences[x] = s;
            }
           
            C_new += s.length();
        }
        cout<<"TIME to stitch: "<<readTimer() - time_a<<" sec."<<endl;
        
        
        
        
        cout<<"GROUP PRINT"<<endl;
        bool* merged = new bool[countNewNode];
        for (int i = 0; i<countNewNode; i++) {
            merged[i] = false;
        }
        
        
        /***MERGE START***/
        if(ALGOMODE == TWOWAYEXT || ALGOMODE == BRACKETCOMP){
            ofstream betterfile;
            betterfile.open("stitchedUnitigs"+modefilename[ALGOMODE]+".fa");
            
            ofstream betterfilePlain;
            betterfilePlain.open("plainOutput"+modefilename[ALGOMODE]+".txt");
            
            
            for ( const auto& p: gmerge.fwdWalkId)
            {
                if(gmerge.fwdVisited[p.first] == false){
                    
                    int fromnode =p.first;
                    int tonode = p.second;
                    deque<int> lst;
                    
                    lst.push_back(fromnode);
                    lst.push_back(tonode);
                    
                    gmerge.fwdVisited[fromnode] = true;
                    gmerge.bwdVisited[tonode] = true;
                    
                    
                    if(gmerge.fwdVisited.count(tonode)>0){
                        while(gmerge.fwdVisited[tonode] == false){
                            gmerge.fwdVisited[tonode] = true;
                            tonode = gmerge.fwdWalkId[tonode];
                            gmerge.bwdVisited[tonode] = true;
                            
                            lst.push_back(tonode);
                            if(gmerge.fwdVisited.count(tonode)==0)
                                break;
                        }
                    }
                    if(gmerge.bwdVisited.count(fromnode)>0){
                        while(gmerge.bwdVisited[fromnode] == false){
                            gmerge.bwdVisited[fromnode] = true;
                            fromnode = gmerge.bwdWalkId[fromnode];
                            gmerge.fwdVisited[fromnode] = true;
                            
                            lst.push_front(fromnode);
                            if(gmerge.bwdVisited.count(fromnode)==0)
                                break;
                        }
                    }
                    
                    string mergeString = "";
                    
                    int headOfThisWalk = walkFirstNode[lst.at(0)]; //CHECK AGAIN
                    assert(!lst.empty());
                    int commonWalkId = lst.at(0);
                    
                    int posOffset = 1;
                    for(auto i: lst){
                        // i is new walk id before merging
                        merged[i] = true;
                        mergeString = plus_strings(mergeString, newSequences[i], K);
                        walkFirstNode[i] = headOfThisWalk;
                        
                        // travesing the walk list of walk ID i
                        for(int uid: newToOld[i]){
                            oldToNew[uid].serial = commonWalkId;
                            oldToNew[uid].finalWalkId = commonWalkId;
                            oldToNew[uid].pos_in_walk = posOffset++;
                            
                        }
                        oldToNew[newToOld[i].back()].isWalkEnd = true;
                        
                    }
                    
                    
                    newNewSequences[headOfThisWalk] = mergeString;
                    
                    
                    //cout<<endl;
                    V_twoway ++;
                    C_twoway+=mergeString.length();
                    betterfile << '>' << commonWalkId <<" LN:i:"<<mergeString.length()<<" ";
                    betterfile<<endl;
                    
                    betterfile<<mergeString;
                    betterfilePlain<<mergeString;
                    
                    betterfile<<endl;
                    betterfilePlain<<endl;
                    
                    
                }
            }
            for (int newNodeNum = 0; newNodeNum<countNewNode; newNodeNum++){
                
                
                if(merged[newNodeNum] == false){
                    oldToNew[newToOld[newNodeNum].back()].isWalkEnd = true;
                    
                    
                    betterfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
                    betterfile<<endl;
                    
                    betterfile<<newSequences[newNodeNum];
                    betterfilePlain<<newSequences[newNodeNum];
                    
                    betterfile<<endl;
                    betterfilePlain<<endl;
                    
                    C_twoway+=newSequences[newNodeNum].length();
                    V_twoway++;
                }
            }
            betterfile.close();
        }
        
        
        /// TWOWAYEXT DONE: NOW LET"S DO BRACK COMP
        
        
        bool* hasSinkconn = new bool[V];
        for (int i = 0; i<V; i++) {
            hasSinkconn[i] = false;
        }
        //@@@@@ BRACKETED
        
        if(1 == 0){
            for (int sinksrc = 0; sinksrc<V; sinksrc++) {
                if(global_issinksource[sinksrc] == 1){
                    list<int> xxx;
                    xxx.push_back(sinksrc);
                    newToOld.push_back(xxx);
                    oldToNew[sinksrc].serial = countNewNode++;
                    oldToNew[sinksrc].finalWalkId = oldToNew[sinksrc].serial;
                    oldToNew[sinksrc].pos_in_walk = 1;
                    oldToNew[sinksrc].isTip = 0;
                    // error resolved in sept 14
                    color[sinksrc] = 'b';
                }
            }
        }
        if(ALGOMODE == BRACKETCOMP){
            if(0==0){
                for (auto const& x : sinkSrcEdges)
                {
                    int sinksrc = x.first;
                    for(edge_t e: x.second){
                        
                        // can this occur?
                        if(color[sinksrc] != 'w'){
                            break;
                        }
                        //there are 3 cases
                        //if consistent this way [[[if(nodeSign[e.toNode] == e.right)]]]
                        //case fwd1: sinksrc -> contig start
                        //case fwd2. sinksrc -> contig middle/end -> ... (say that sinksrc is LEFT)
                        //case fwd3. sinksrc -> sinksrc_other (i'd say ignore this for now)
                        //
                        
                        //case bwd1. contig end -> sinksrc
                        //case bwd2. .... -> contig middle/start -> sinksrc (say that sinksrc is RIGHT)
                        //case bwd3. sinksrc_other -> sinksrc  (i'd say ignore this for now)
                        
                        // 3 fwd cases
                        if(nodeSign[e.toNode] == e.right){  //ensure it is a fwd case
                            if(color[e.toNode]!='w' && color[e.toNode]!='r' && color[e.toNode]!='l'){//  this ensures that this to vertex is NOT sinksrc_other
                                //case 1, 2
                                int whichwalk = oldToNew[e.toNode].finalWalkId;
                                //*** case fwd1 : sinksrc -> contigStart
                                //case fwd2. sinksrc -> contig middle/end -> ... (say that sinksrc is LEFT)
                                //let's merge case fwd1 & fwd2
                                //color[sinksrc] = 'b';
                                
                                nodeSign[sinksrc] = e.left;
                                color[sinksrc] = 'l';
                                oldToNew[sinksrc].serial = whichwalk;
                                oldToNew[sinksrc].finalWalkId = whichwalk;
                                oldToNew[sinksrc].pos_in_walk = oldToNew[e.toNode].pos_in_walk;
                                
                                //oldToNew[sinksrc].isTip = -1; // from left
                                
                                if(oldToNew[e.toNode].pos_in_walk == 1 && !hasSinkconn[e.toNode] ){
                                    oldToNew[sinksrc].isTip = 0;
                                    hasSinkconn[e.toNode] = true;
                                    oldToNew[sinksrc].pos_in_walk--;
                                }else{
                                    oldToNew[sinksrc].isTip = -1;
                                }
                                
                                
                                
                                
                                //fwd1
                                //newToOld[whichwalk].insert(newToOld[whichwalk].begin(), sinksrc);
                                //fwd2
                                //std::list<int>::iterator it;
                                //it = find (newToOld[whichwalk].begin(), newToOld[whichwalk].end(), e.toNode);
                                //newToOld[whichwalk].insert(it, sinksrc);
                                
                            }
                            
                        }else{
                            // 3 bwd cases
                            
                            if(color[e.toNode]!='w' && color[e.toNode]!='r' && color[e.toNode]!='l'){
                                int whichwalk = oldToNew[e.toNode].finalWalkId;
                                
                                //*** case bwd1: contigend --> sinksrc
                                //*** case bwd2: contigmiddle--> sinksrc
                                
                                
                                nodeSign[sinksrc] = !e.left;
                                //color[sinksrc] = 'b';
                                color[sinksrc] = 'r';
                                oldToNew[sinksrc].serial = whichwalk;
                                oldToNew[sinksrc].finalWalkId = whichwalk;
                                oldToNew[sinksrc].pos_in_walk = oldToNew[e.toNode].pos_in_walk ;
                                
                                
                                if(oldToNew[e.toNode].isWalkEnd == true && !hasSinkconn[e.toNode] ){
                                    oldToNew[sinksrc].isTip = 0;
                                    hasSinkconn[e.toNode] = true;
                                    oldToNew[sinksrc].pos_in_walk++;
                                }else{
                                    oldToNew[sinksrc].isTip = 1;
                                }
                                //oldToNew[sinksrc].isTip = 1; //from right
                                
                                //C_new+= unitigs.at(sinksrc).ln - (K-1) + 2;
                                //C_bracketed+= unitigs.at(sinksrc).ln - (K-1) + 2;
                                //newToOld[whichwalk].insert(newToOld[whichwalk].end(), sinksrc);
                                
                            }
                        }
                    }
                }
            }
            
            // now take care of all the remaining edges
//            for (auto const& x : sinkSrcEdges)
//            {
//                int sinksrc = x.first;
//                if(color[sinksrc] == 'w'){  //still white, that means it goes isolated now
//                    list<int> xxx;
//                    xxx.push_back(sinksrc);
//                    newToOld.push_back(xxx);
//                    oldToNew[sinksrc].serial = countNewNode++;
//                    oldToNew[sinksrc].finalWalkId = oldToNew[sinksrc].serial;
//                    oldToNew[sinksrc].pos_in_walk = 1;
//                    oldToNew[sinksrc].isTip = 0;
//                    // error resolved in sept 14
//                    color[sinksrc] = 'b';
//                }
//            }
            
            for (int sinksrc = 0; sinksrc<V; sinksrc++) {
                if(global_issinksource[sinksrc] == 1 && color[sinksrc] == 'w' ){
                    list<int> xxx;
                    xxx.push_back(sinksrc);
                    newToOld.push_back(xxx);
                    oldToNew[sinksrc].serial = countNewNode++;
                    oldToNew[sinksrc].finalWalkId = oldToNew[sinksrc].serial;
                    oldToNew[sinksrc].pos_in_walk = 1;
                    oldToNew[sinksrc].isTip = 0;
                    // error resolved in sept 14
                    color[sinksrc] = 'b';
                }
            }
            
            
            
            
            //BRACKETCOMP encoder and printer::::
            vector<mytuple> sorter;
            for(int uid = 0 ; uid< V; uid++){
                new_node_info_t nd = oldToNew[uid];
                //if(!global_issinksource[uid]){
                    sorter.push_back(make_tuple(uid, nd.finalWalkId, nd.pos_in_walk, nd.isTip));
                //}
                
            }
            stable_sort(sorter.begin(),sorter.end(),sort_by_tipstatus);
            stable_sort(sorter.begin(),sorter.end(),sort_by_pos);
            stable_sort(sorter.begin(),sorter.end(),sort_by_walkId);
            
            /// START OUTPUTTING
            ofstream tipFile;
            tipFile.open("tipOutput.txt");
            
            ofstream tipDebugFile;
            tipDebugFile.open("tipDebug.txt");
            
            
            int lastWalk = -1;
            string walkString = "";
            string tipLessWalkString ="";
            
            for(mytuple n : sorter){
                int uid = get<0>(n);
                int finalWalkId = get<1>(n);
                int pos_in_walk = get<2>(n);
                int isTip = get<3>(n);
                //cout<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;

                string unitigString;
                if(finalWalkId!=lastWalk){
                    if(lastWalk != -1){
                        //print previous walk
                        tipDebugFile<<">"<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
                        tipFile<< '>' << lastWalk << "\n" ;
                        V_bracketed++;
                        C_bracketed+=walkString.length();
                        
                        tipDebugFile<<walkString<<endl;
                        tipFile<< walkString<<endl;
                        //cout<<endl;
                    }
                    
                    //start a new walk
                   // cout<<"Walk: (" <<finalWalkId<<" ) = ";
                    walkString = "";
                    lastWalk = finalWalkId;
                    
                }else{
                    //cout<<pos_in_walk<<" ";
                }
                
                if(nodeSign[uid] == false){
                    unitigString =  reverseComplement(unitigs.at(uid).sequence);
                }else{
                    unitigString =  (unitigs.at(uid).sequence);
                }
                
               
                if(isTip == 0){
                     walkString = plus_strings(walkString, unitigString, K);
                }else if(isTip==1){ //right R   R    ]   ]   ]   ]
                    //cut prefix
                    if(0==0){
                        unitigString = unitigString.substr(K - 1, unitigString.length() - (K - 1));
                        if(walkString.length()<K){
                            cout<<"pos: "<<walkString.length()<<endl;
                        }
                        walkString += "]" + unitigString + "]";
                    }
                    if(0==1){
                        tipFile<<">pref\n"<<unitigString<<endl;
                    }
                    
                }else if(isTip==-1){ //left L   L    [ [ [
                    //cut suffix
                    if(0==0){
                        unitigString = unitigString.substr(0, unitigString.length() - (K - 1));
                        if(walkString.length()<K){
                            cout<<"pos: "<<walkString.length()<<endl;
                        }
                        walkString += "[" + unitigString + "[";
                        
                    }
                    if(1==0){
                         tipFile<<">suf\n"<<unitigString<<endl;
                    }
                    
                }
                    
                    
                tipDebugFile<<">"<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
                tipFile<<">"<<lastWalk<<endl;
                
            }
            V_bracketed++;
            C_bracketed+=walkString.length();
            
            tipDebugFile<< walkString;
            tipDebugFile<<endl;
            tipDebugFile.close();
            
            
            tipFile<< walkString;
            tipFile<<endl;
            tipFile.close();
            
            //DECODER
        }
        
        
        
        
        
        
        delete[] merged;
        delete []  global_issinksource;
        delete []  global_priority;

    }

    ~Graph() {
        delete [] color;
        delete [] p;
        delete [] nodeSign;
        delete [] oldToNew;
        delete [] saturated;
        delete [] indegree;
        delete [] countedForLowerBound;
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


void formattedOutputForwardExt(Graph &G){
    string plainOutput = "plainOutput"+modefilename[ALGOMODE]+".txt";
    ofstream plainfile;
    plainfile.open(plainOutput);
    
    string stitchedUnitigs = "stitchedUnitigs"+modefilename[ALGOMODE]+".fa";
    ofstream myfile;
    myfile.open (stitchedUnitigs);
    //>0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:-
    for (int newNodeNum = 0; newNodeNum<G.countNewNode; newNodeNum++){
        myfile << '>' << newNodeNum <<" LN:i:"<<newSequences[newNodeNum].length()<<" ";
        //plainfile << '>' << newNodeNum;
        //C_new+=newSequences[newNodeNum].length();
        //plainfile<<endl;
        myfile<<endl;
        
        plainfile<<newSequences[newNodeNum];
        myfile<<newSequences[newNodeNum];
        
        plainfile<<endl;
        myfile<<endl;
    }
    //myfile << '>' << newNodeNum <<">0 LN:i:13 KC:i:12 km:f:1.3  L:-:0:- L:-:2:-  L:+:0:+ L:+:1:- " ;
    myfile.close();
    plainfile.close();
    
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

int get_data(const string& unitigFileName,
             uchar*& data,
             vector<unitig_struct_t>& unitigs,
             uint64_t& char_count
             ) {
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
    
    cout << "Complete reading input." << endl;
    
    return EXIT_SUCCESS;
}

pair<int, int> getFileSizeKB(string fn = "plainOutput.txt"){
    pair<int, int> sizeBeforeAfter;
    
    string fngz = fn+".gz"; //"plainOutput.fa.gz";
    string countFile = "incount.txt";
    ostringstream stringStream;
    string copyOfStr;
    string line;
    
    system(("du -k " + fn + " | cut -f1 > "  + countFile).c_str()); // du -h plainOutput.fa | cut -f1 > incount.txt
    ifstream cf;
    cf.open(countFile);
    getline(cf, line);
    line = delSpaces(line);
    sscanf(line.c_str(), "%d", &sizeBeforeAfter.first);
    cf.close();
    
  
    
    system(("gzip "+fn + " \n" + " du -k " + fngz + " | cut -f1 > "  + countFile).c_str()); // du -h plainOutput.fa.gzip | cut -f1 > incount.txt
    cf.open(countFile);
    getline(cf, line);
    line = delSpaces(line);
    sscanf(line.c_str(), "%d", &sizeBeforeAfter.second);
    cf.close();
    
    //sizeBeforeAfter.first *= 1024*8;
    //sizeBeforeAfter.second *= 1024*8;
    system("rm -rf incount.txt");
    return sizeBeforeAfter;
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

int main(int argc, char** argv) {
    FILE * statFile;
    statFile = fopen (("stats"+modefilename[ALGOMODE]+".txt").c_str(),"w");

//    string debugFileName = "debug.txt";
//    ofstream debugFile;
//    debugFile.open(debugFileName);
    
    
    const char* nvalue = "" ;
    
    int c ;
    
    ///*
    while( ( c = getopt (argc, argv, "i:k:m:d:f:") ) != -1 )
    {
        switch(c)
        {
            case 'i':
                if(optarg) nvalue = optarg;
                break;
            case 'f':
                if(optarg) {
                    FLG_NEWUB = static_cast<bool>(std::atoi(optarg));
                }
                break;
            case 'm':
                if(optarg) {
                    ALGOMODE = static_cast<ALGOMODE_T>(std::atoi(optarg));
                }
                break;
            case 'd':
                if(optarg) {
                    DBGFLAG = static_cast<DEBUGFLAG_T>(std::atoi(optarg));
                }
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
    
    if(K==0 || strcmp(nvalue, "")==0){
        fprintf(stderr, "Usage: %s -k <kmer size> -i <input-file-name>\n",
                argv[0]);
        exit(EXIT_FAILURE);
    }
    
    
    
    UNITIG_FILE = string(nvalue);
    //*/
    
    ifstream infile(UNITIG_FILE);
    if(!infile.good()){
        fprintf(stderr, "Error: File named \"%s\" cannot be opened.\n", UNITIG_FILE.c_str());
        exit(EXIT_FAILURE);
    }
    
    
    uint64_t char_count;
    uchar *data = NULL;
    
    double startTime = readTimer();
    cout << "## START reading file: " << UNITIG_FILE << ": K = "<<K<<endl;
    if (EXIT_FAILURE == get_data(UNITIG_FILE, data, unitigs, char_count)) {
        return EXIT_FAILURE;
    }
    infile.close();
    double TIME_READ_SEC = readTimer() - startTime;
    cout<<"TIME to read file "<<TIME_READ_SEC<<" sec."<<endl;
    
    
    Graph G;
    
    //count total number of edges
    int E = 0;
    for (int i = 0; i < G.V; i++) {
        E += adjList[i].size();
    }
    int V = G.V;
    int numKmers = 0;
    int C = 0;
    
    for (unitig_struct_t unitig : unitigs) {
        C += unitig.ln;
        numKmers +=  unitig.ln - K + 1;
    }
    
//    if(DBGFLAG == NODENUMBER_DBG){
//        cout<<"Total Nodes: "<<V<<" Edges: "<<E<<" K-mers: "<<numKmers<<endl;
//    }
    
    cout<<"## START gathering info about upper bound. "<<endl;
    double time_a = readTimer();
    G.indegreePopulate();
    
    cout<<"TIME for information gather: "<<readTimer() - time_a<<" sec."<<endl;
    
    delete [] global_indegree;
    delete [] global_outdegree;
    delete [] global_plusindegree;
    delete [] global_plusoutdegree;
    
    
    if(ALGOMODE == GRAPHPRINT){
        char sss[1000];
        cout<<"input the nodes to include in printing separated by space (i.e. 20 19 18): "<<endl;
        while(true){
            //string ipstr = "20 19 18";
            gets(sss);
            string ipstr(sss);
            if(ipstr=="stop"){
                break;
            }
            makeGraphDot(ipstr);
            cout<<"done print, say again:"<<endl;
        }
    }
    
    
    int walkstarting_node_count = ceil((sharedparent_count + sink_count + source_count)/2.0) + isolated_node_count;
    int charLowerbound = C-(K-1)*(G.V - walkstarting_node_count*1.0);
    float upperbound = (1-((C-(K-1)*(G.V - walkstarting_node_count*1.0))/C))*100.0;
    
    printf( "%d\t\
           %d\t\
           %d\t\
           %d\t\
           %d\t\
           %d\t\
           %.2f\t\
           %.2f%%\t\
           %d\t\
           %d\t\
           %d\t\
           %d\t\
           %.2f%%\t",
           K,
           numKmers,
           V,
           E,
           C,
           charLowerbound,
           (charLowerbound*2.0)/numKmers,
           upperbound,
           isolated_node_count,
           sink_count,
           source_count,
           sharedparent_count,
           sharedparent_count*100.0/V
           );
    
    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
        printf("%.2f%%\t", (i->second)*100.0/V);
    }
    printf("%.2f%%\t\
           %.2f%%\t",
           isolated_node_count*100.0/V,
           (sink_count+source_count)*100.0/V);
    //fprintf(statFile, "\n");

    // Iterating the map and printing ordered values
//    for (auto i = inOutCombo.begin(); i != inOutCombo.end(); i++) {
//        cout << "(" << i->first.first<< ", "<< i->first.second << ")" << " := " << i->second << '\n';
//    }
    
    if(ALGOMODE == PROFILE_ONLY){
        printf("\n");
        return 0;
    }

    
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
//##################################################//
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@//
//##################################################//
    
    //STARTING DFS
    cout<<"## START DFS: "<<endl;
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
 
    
    //Stats after all done
    //int V_new = G.countNewNode;
    //int C_new = 0;
    
    double TIME_TOTAL_SEC = readTimer() - startTime;
//    map<int, string>::iterator it;
//    int maxlen = 0;
//    for (it = newSequences.begin(); it != newSequences.end(); it++)
//    {
//        C_new += (it->second).length();
//        if((it->second).length() >maxlen){
//            maxlen =(it->second).length();
//        }
//    }
    

    // For collecting stats
    //int U_MAX = maximumUnitigLength();
    
    
    
    time_a = readTimer();
    
    
    if(ALGOMODE==BRACKETCOMP){
        //formattedOutputForwardExt(G);
        C_new = C_bracketed;
        V_new  = V_bracketed;
        //V_new = G.countNewNode;
    }else if(ALGOMODE == TWOWAYEXT){
        V_new = V_twoway;
        C_new = C_twoway;
    }else if(ALGOMODE == BASIC){
        formattedOutputForwardExt(G);
        V_new = G.countNewNode;
        //C_new = C_new;
    }else{
        formattedOutputForwardExt(G);
        V_new = G.countNewNode;
        //C_new = C_new;
    }
    
    cout<<"TIME to output: "<<readTimer() - time_a<<" sec."<<endl;
    
    
    
    
    
   float percent_saved_c = (1-(C_new*1.0/C))*100.0;
    float theoreticalBitsKmerSaved;
    if(ALGOMODE== BRACKETCOMP){
         theoreticalBitsKmerSaved =C_new*3.0/numKmers;
    }else{
         theoreticalBitsKmerSaved =C_new*2.0/numKmers;
    }
   
    
    printf("%s\t",  mapmode[ALGOMODE].c_str());
    printf("%d\t\
            %.2f%%\t\
            %.2f%%\t\
           %d\t\
           %.2f\t",
           V_new,
           percent_saved_c,
            upperbound - percent_saved_c,
           C_new,
           theoreticalBitsKmerSaved
            );
    printf("%.2f\t\
           %.2f\t",
           TIME_READ_SEC,
           TIME_TOTAL_SEC
           );
//    fprintf(statFile, "%.2f\t\
//           %.2f\t",
//           getFileSizeKB().first*1.0/numKmers*1024*8,
//           getFileSizeKB().second*1.0/numKmers*1024*8);
    printf("\n");
    printf("\n");
    
    fclose(statFile);
    return EXIT_SUCCESS;
}
