//
//  absorbGraph.hpp

#ifndef absorbGraph_hpp
#define absorbGraph_hpp

#include "global.hpp"

//make the absorb graph with cycle
void sorterIndexAndAbsorbGraphMaker(vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap, map<int, vector<edge_t> > & absorbGraph, char*& absorbedCategory){
    int prevWalkId = -1;
    int lastWalkStartingIndex = -1;
    
    for(int tup_i = 0; tup_i < sorter.size(); tup_i++){
        
        MyTypes::fourtuple tup = sorter[tup_i];
        int uid = get<0>(tup);
        new_node_info_t nd = oldToNew[uid];
        int finalWalkId = get<1>(tup);
        
        // ideally we should only keep one of them, either add walk start or walk end
        if(prevWalkId !=finalWalkId ){  //walk starting
            lastWalkStartingIndex = tup_i;
            sorterIndexMap[finalWalkId] = lastWalkStartingIndex;
            
            vector<edge_t> adju = adjList.at(uid);
            for (edge_t e : adju) {
                int absorberWalk = oldToNew[e.toNode].finalWalkId;
                int absorbedWalk = nd.finalWalkId;

                if(absorberWalk != absorbedWalk && oldToNew[e.toNode].pos_in_walk!=1){
                    //the edge "enew" is from absorber uid to absorbed uid
                    edge_t enew;
                    if(e.left==e.right){
                        enew.left = (e.left?false:true);
                        enew.right = (e.right?false:true);
                    }else{
                        enew.left = e.left;
                        enew.right = e.right;
                    }
                    
                    enew.toNode = uid;
                    absorbGraph[e.toNode].push_back(enew);
                    
                    //cout<<absorberWalk<<"->absorbs->"<<absorbedWalk<<endl;
                    //TODO: need to fix for decoding
                    //cout<<"hey"<<endl;
                    //break;
                }
            }
        }
        
        if(1==0){
            if(tup_i + 1 < sorter.size()){
                if(get<1>(sorter[tup_i + 1])!=finalWalkId && lastWalkStartingIndex!=tup_i){   //so it is the end of a walk, and it is not an isolated walk (with one vertex only)
                    vector<edge_t> adju = adjList.at(uid);
                    for (edge_t e : adju) {
                        int absorberWalk = oldToNew[e.toNode].finalWalkId;
                        int absorbedWalk = nd.finalWalkId;
                        if(absorberWalk != absorbedWalk){
                            //absorbGraph[e.toNode].push_back(enew);
                        }
                    }
                }
            }
        }
        prevWalkId = finalWalkId;
    }
}



void absorbGraphIsCyclicUtil(int uid, char*& visited, map<int, vector<edge_t> >& absorbGraph, map<int, stack<edge_t> >& absorbGraphCycleRemoved, queue<int>& orderOfUnitigs, char*& absorbedCategory )
{
    /**
     map<int, int > sorterIndexMap: map walkID to location of <walk start, end> in sorter
     //v = visited[oldToNew[uid].finalWalkId] is the absorber
     ****/

    visited[oldToNew[uid].finalWalkId] = 'g';
    vector<edge_t> adjv = absorbGraph[uid];
    
    for(edge_t e : adjv){
        int i = e.toNode;
        if (visited[oldToNew[i].finalWalkId] == 'g'){
            //removed edge
        }else if (visited[oldToNew[i].finalWalkId] == 'w'){
            //add to the list                //find absorption category
            //call the absorber v[i]
            //.....orderOfUnitigs.push(i);
            
            absorbGraphCycleRemoved[uid].push(e);
            
            if(e.left == nodeSign[uid] && e.right == nodeSign[e.toNode] )
                absorbedCategory[e.toNode] = '1';
            if(e.left == nodeSign[uid] && e.right != nodeSign[e.toNode] )
                absorbedCategory[e.toNode] = '2';
            if(e.left != nodeSign[uid] && e.right != nodeSign[e.toNode] )
                               absorbedCategory[e.toNode] = '3';
            if(e.left != nodeSign[uid] && e.right == nodeSign[e.toNode] )
                                   absorbedCategory[e.toNode] = '4';
            absorbGraphIsCyclicUtil(i, visited, absorbGraph, absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
        }
    }
    visited[oldToNew[uid].finalWalkId]  = 'b';
    return;
}

void  removeCycleFromAbsorbGraph(vector<MyTypes::fourtuple>& sorter,  map<int, int >& sorterIndexMap, map<int, vector<edge_t> >& absorbGraph, map<int, stack<edge_t> >& absorbGraphCycleRemoved,  queue<int>& orderOfUnitigss, char*& absorbedCategory )
{
   // int Vuid = absorbGraph.size();
   // int Vwalkid = sorterIndexMap.size();  // only the number of walks participating
    char *visited = new char[countNewNode];
    
    //
    //for(int i = 0; i < Vwalkid; i++)
    for(int i = 0; i < countNewNode; i++)
    {
        if(obsoleteWalkId[i]){
            visited[i] = 'b';
        }else{
             visited[i] = 'w';
        }
       
    }
  
    // Call the recursive helper function to detect cycle in different DFS trees
    for(int i = 0; i <adjList.size(); i++){
        if(absorbGraph.count(i)>0){
            //for(int i = 0; i <adjList.size(); i++){
            if(visited[oldToNew[i].finalWalkId] == 'w'){
                orderOfUnitigss.push(i);
                absorbGraphIsCyclicUtil(i, visited, absorbGraph, absorbGraphCycleRemoved, orderOfUnitigss, absorbedCategory);
            }
        }
    
    }
    //map int, priority queue uid loc
}

string recursiveWalkStringMaker(int& startWalkIndex, vector<bool>& isItAPrintedWalk, vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap,  map<int, stack<edge_t> >& absorbGraphCycleRemoved, queue<int>& orderOfUnitigs, char*& absorbedCategory){
    
    bool isThisAbsorbedWalk = false;;
    string unitigString;
    string walkString = "";
    
    while(true){
        assert(startWalkIndex<sorter.size());
        MyTypes::fourtuple n = sorter[startWalkIndex];
       int finalWalkId = get<1>(n);
       
       //@ABSORB
       if (isItAPrintedWalk[finalWalkId]) return "";

       int uid = get<0>(n);
        
//        if(uid==4869){
//            cout<<"";
//        }
        
//        if(finalWalkId==11){
//            cout<<"";
//        }
        
       int pos_in_walk = get<2>(n);
       int isTip = get<3>(n);
       //cout<<uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
        
        if(nodeSign[uid] == false){
            unitigString =  reverseComplement(unitigs.at(uid).sequence);
        }else{
            unitigString =  (unitigs.at(uid).sequence);
        }
        
        stack<edge_t> stType12;
        stack<edge_t> stType34;
        if(absorbGraphCycleRemoved.count(uid) > 0){         /*populate two types of stacks*/
            stack<edge_t> st = absorbGraphCycleRemoved[uid];
            while(!st.empty()){
                edge_t st_top = st.top();
                if(absorbedCategory[st_top.toNode]=='3' or absorbedCategory[st_top.toNode]=='4'  ){
                    stType34.push(st_top);
                }else if(absorbedCategory[st_top.toNode]=='1' or absorbedCategory[st_top.toNode]=='2' ){
                    stType12.push(st_top);
                }else{
                    cout<<"errorrrrrrrrrrr"<<endl;
                    assert(false);
                }
                st.pop();
            }
        }
        
        if(pos_in_walk == 1 && absorbedCategory[uid]=='0'){
            walkString+= cutSuf(unitigString, K);

            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                stType34.pop();
                walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
            }
            
            walkString+= suf(unitigString, K);
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                stType12.pop();
                walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
            }
        }
        
        
        if(pos_in_walk != 1 && absorbedCategory[uid]=='0'){
            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                stType34.pop();
                walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
            }
            
            walkString+= cutPref(unitigString, K);
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                stType12.pop();
                walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
            }
        }
        if(pos_in_walk != 1 && absorbedCategory[uid]!='0'){
            assert(false);
        }
        
        if(pos_in_walk == 1 && absorbedCategory[uid]!='0'){
            isThisAbsorbedWalk=true;
            if(absorbedCategory[uid]=='2' or absorbedCategory[uid]=='3'){
                walkString+= cutSuf(unitigString, K);
            }
            if(absorbedCategory[uid]=='1') walkString+="+";
            else if(absorbedCategory[uid]=='4') walkString+="-";
            while(!stType34.empty()){
                int sindex = sorterIndexMap[oldToNew[stType34.top().toNode].finalWalkId];
                stType34.pop();
                walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
            }
            
            if(absorbedCategory[uid]=='1' or absorbedCategory[uid]=='4'){
                           walkString+= cutPref(unitigString, K);
                       }
            
            while(!stType12.empty()){
                int sindex = sorterIndexMap[oldToNew[stType12.top().toNode].finalWalkId];
                stType12.pop();
                walkString+= recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
            }
            if(absorbedCategory[uid]=='3') walkString+="+";
            else if(absorbedCategory[uid]=='2') walkString+="-";
        }


    
        if(1==0 && isTip == 0){//
        
            if(absorbGraphCycleRemoved.count(uid) > 0){
                stack<edge_t> st = absorbGraphCycleRemoved[uid];
                while(!st.empty()){
                    int st_top_uid = st.top().toNode;
                    if(1==1){
                        int sindex = sorterIndexMap[oldToNew[st_top_uid].finalWalkId];
                        st.pop();

                        string absorbedWalk = recursiveWalkStringMaker(sindex, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
                        //depending on the type of edge
                        if(absorbedWalk!=""){
                            if(absorbedCategory[st_top_uid] == '1'){
                                walkString += "{+" + absorbedWalk.substr(K - 1, absorbedWalk.length() - (K - 1)) + "}";
                            }else if(absorbedCategory[st_top_uid] == '2'){
                                walkString += "{-" + absorbedWalk.substr(0, unitigString.length() - (K - 1)) + "}";
                                
                            }else if(absorbedCategory[st_top_uid] == '3'){
                                walkString += "(+" + absorbedWalk.substr(K - 1, absorbedWalk.length() - (K - 1)) + ")";
                            }else if(absorbedCategory[st_top_uid] == '4'){
                                walkString += "(-" + absorbedWalk.substr(K - 1, absorbedWalk.length() - (K - 1)) + ")";
                            }
                        }
                    }else{
                        //st.pop();
                    }
                }
            }
        }else if(isTip==1 && MODE_ABSORPTION_TIP){ //right R   R    ]   ]   ]   ]
            //cut prefix: correct
            if(0==0){
                unitigString = unitigString.substr(K - 1, unitigString.length() - (K - 1));
                if(walkString.length()<K){
                    cout<<"pos: "<<walkString.length()<<endl;
                }
                walkString += "(" + unitigString + ")";
            }
            if(1==0){
                tipFile<<">pref\n"<<unitigString<<endl;
            }
            
        }else if(isTip==2 && MODE_ABSORPTION_TIP){ //left L   L    [ [ [
            //cut suffix: correct
            if(0==0){
                unitigString = unitigString.substr(0, unitigString.length() - (K - 1));
                if(walkString.length()<K){
                    cout<<"pos: "<<walkString.length()<<endl;
                }
                walkString += "{" + unitigString + "}";
            }
            if(1==0){
                tipFile<<">suf\n"<<unitigString<<endl;
            }
        }
            
        if(startWalkIndex+1 == sorter.size()) {
            //print previous walk
            tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
            // V_tip_ustitch++;
           // C_tip_ustitch+=walkString.length();
//
            tipDebugFile<<walkString<<endl;
          //  tipFile<< walkString<<endl;
            
            
            isItAPrintedWalk[finalWalkId] = true;
            break;
        }else if(get<1>(sorter[startWalkIndex+1]) != finalWalkId){
                        //print previous walk
                        tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
            //           V_tip_ustitch++;
            //            C_tip_ustitch+=walkString.length();
            //
                        tipDebugFile<<walkString<<endl;
             //           tipFile<< walkString<<endl;
                        
                        isItAPrintedWalk[finalWalkId] = true;
                        break;
        }
        startWalkIndex++;
    }
    if(isThisAbsorbedWalk) walkString="["+walkString+"]";
    return walkString;
}

void tipAbsorbedOutputter(vector<MyTypes::fourtuple>& sorter, map<int, int >& sorterIndexMap, map<int, stack<edge_t> >& absorbGraphCycleRemoved, queue<int>& orderOfUnitigs, char*& absorbedCategory ){
    /// START OUTPUTTING
    vector<bool> isItAPrintedWalk;
    for(int i = 0; i< countNewNode; i++) isItAPrintedWalk.push_back(false);

    tipFile.open("tipOutput.txt");
    tipDebugFile.open("tipDebug.txt");
    
    while(!orderOfUnitigs.empty()){
        int it = sorterIndexMap[oldToNew[orderOfUnitigs.front()].finalWalkId];
        
        if(DBGFLAG==PRINTER){
            cout<<"order: "<<orderOfUnitigs.front()<<endl;
        }
        
        orderOfUnitigs.pop();
        string walkString = recursiveWalkStringMaker(it, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
        //tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
        
        if(walkString!=""){
            V_tip_ustitch++;
            C_tip_ustitch+=walkString.length();
            
           // tipDebugFile<<walkString<<endl;
            tipFile<< walkString<<endl;
        }
    }
    
    //exit(1);
    
    int it = -1;
    while(true){
        if(orderOfUnitigs.empty()){
            it++;
            if(it>=sorter.size()) break;
        }else{
            //n = sorter[sorterIndexMap[orderOfUnitigs.front()]];
            it = sorterIndexMap[oldToNew[orderOfUnitigs.front()].finalWalkId];
            orderOfUnitigs.pop();
        }
        //recursive call
        string walkString = recursiveWalkStringMaker(it, isItAPrintedWalk, sorter, sorterIndexMap,  absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
        //tipDebugFile<<">"<<finalWalkId << " " << uid<<" " <<finalWalkId<<" "<<pos_in_walk<<" "<<isTip<<endl;
        
        if(walkString!=""){
            V_tip_ustitch++;
            C_tip_ustitch+=walkString.length();
            
           // tipDebugFile<<walkString<<endl;
            tipFile<< walkString<<endl;
        }
        //
    }
    //delete [] isItAnAbsorbedWalk;
    tipDebugFile.close();
    tipFile.close();
}

//map<int, vector<edge_t> >
void print_absorb_graph(map<int,vector<edge_t> > const &m)
{
    cout<<"ABSORB GRAPH---------------"<<endl;
    for (auto const& pair: m) {
        std::cout << pair.first << ":  ";
        vector<edge_t> vec = pair.second;
        for(int i=0; i<vec.size(); i++){
            edge_t e = vec[i];
            cout<< e.toNode  << ", ";
        }
        cout<<"\n";
    }
    cout<<"ABSORB GRAPH END---------------"<<endl;
}

//map<int, vector<edge_t> >
void print_absorb_graph_acyclic(map<int,stack<edge_t> > const &m)
{
    cout<<"ABSORB GRAPH ACYCLIC---------------"<<endl;
    for (auto const& pair: m) {
        std::cout << pair.first << ":  ";
        stack<edge_t> vec = pair.second;
        while(!vec.empty()){
            edge_t e = vec.top();
            vec.pop();
            cout<< e.toNode  << ", ";
        }
        cout<<"\n";
    }
    cout<<"ABSORB GRAPH ACYCLIC END---------------"<<endl;
}



void absorptionManager(vector<MyTypes::fourtuple> sorter) {

    if(DBGFLAG==PRINTER){
        //    int uid = get<0>(n);
        //    int finalWalkId = get<1>(n);
        //    int pos_in_walk = get<2>(n);
        for(auto i: sorter ){
            cout<<"walk "<<get<1>(i)<<":";
            cout<<"("<<get<2>(i)<<")";
             cout<<get<0>(i)<<"";
            cout<<endl;
        }
    }
    
    int totuids = adjList.size();
    char* absorbedCategory = new char[adjList.size()];
    for(int i=0; i< totuids; i++){
        absorbedCategory[i] = '0';
    }
    
    map<int, int > sorterIndexMap;
    map<int, vector<edge_t> > absorbGraph;
    if(MODE_ABSORPTION_NOTIP){
        sorterIndexAndAbsorbGraphMaker(sorter, sorterIndexMap, absorbGraph, absorbedCategory);
    }
    
   
    if(DBGFLAG==PRINTER){
        print_absorb_graph(absorbGraph);
    }
    
    
    map<int, stack<edge_t> > absorbGraphCycleRemoved;
    

    queue<int> orderOfUnitigs;
     if(MODE_ABSORPTION_NOTIP ){
    removeCycleFromAbsorbGraph(sorter,  sorterIndexMap, absorbGraph, absorbGraphCycleRemoved,  orderOfUnitigs,  absorbedCategory);//last two are output
    absorbGraph.clear();
     }
    
    if(DBGFLAG==PRINTER){
        for (int i=0; i<totuids; i++) {
            cout<<i<<"->"<<absorbedCategory[i]<<endl;
        }
         print_absorb_graph_acyclic(absorbGraphCycleRemoved);
    }
    
    tipAbsorbedOutputter(sorter, sorterIndexMap, absorbGraphCycleRemoved, orderOfUnitigs, absorbedCategory);
    
    
    delete [] absorbedCategory;
}

#endif /* absorbGraph_hpp */
