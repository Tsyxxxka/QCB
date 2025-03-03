#include <bits/stdc++.h>
#define PII pair<int,int>
using namespace std;
int s,t;
unordered_map<int, PII> id2edge;
unordered_map<int, unordered_map<int, int>> edge2id;
vector<int> node_scc;
vector<unordered_set<int>> cycle_basis;
vector<int> start_cb;
vector<vector<int>> sup_graph,edge_cycle;
vector<unordered_set<int>> join_graph;
int n, m, cb_num;
int edge_st;
vector<int> low, dfn, in_stack, ss, scc, sz;
int dfncnt = 1, scc_num = 0, tp = 0;
int circuit_cnt = 0;
unordered_set<int> biedge;

vector<int> join_visited, start_cycle;
unordered_set<int> current_cir, connect_cycle;
double t2;
int connect_component=0;
unordered_set<int> all_nodes;
int ttt=1;

void join_cycles(int new_cid){
    for (auto new_eid: cycle_basis[new_cid]){
        if (current_cir.find(new_eid)!=current_cir.end()) current_cir.erase(new_eid);
        else current_cir.insert(new_eid);
    }
}

unordered_map<int,unordered_set<int>> judge_cycle;
unordered_map<int,int> cycle_visited;

void dfs_cycle(int root){
    cycle_visited[root] = connect_component;
    for (auto rn: judge_cycle[root]) {
        if (cycle_visited[rn]==0) {
            dfs_cycle(rn);
        }
    }
}

bool judge(){
    connect_component=0;
    all_nodes.clear();
    judge_cycle.clear(); 
    int f1,f2;
    for (auto eid: current_cir) {
        f1 = id2edge[eid].first;
        f2 = id2edge[eid].second;
        all_nodes.insert(f1); all_nodes.insert(f2);
        unordered_set<int> new_set; new_set.clear();
        if (judge_cycle.find(f1)!=judge_cycle.end()) judge_cycle[f1].insert(f2);
        else judge_cycle[f1] = new_set, judge_cycle[f1].insert(f2);
        if (judge_cycle.find(f2)!=judge_cycle.end()) judge_cycle[f2].insert(f1);
        else judge_cycle[f2] = new_set, judge_cycle[f2].insert(f1); 
    }
    cycle_visited.clear();
    for (auto kk: all_nodes) {
        cycle_visited[kk]=0;
    }
    for (auto nn: all_nodes) {
        if (cycle_visited[nn]==0) {
            connect_component++;
            dfs_cycle(nn);
        }
    }
    if (connect_component>1) {
        return false;
    }
    return true;
}

vector<unordered_set<int>> cycle_component_cb_set; 

void search_full_cycle(int is_new, unordered_set<int> current_nb){
    if (is_new) {
        if (judge()) {
            circuit_cnt++;
            if (circuit_cnt%100==0) {
                double t3=clock();
                cout << circuit_cnt << "\t" << (t3 - t2) / CLOCKS_PER_SEC << endl;
            }
        }
        else {
            assert(connect_component>1);
            cycle_component_cb_set.clear();
            cycle_component_cb_set.resize(connect_component);
            for (auto n_e: judge_cycle) {
                for (auto se: n_e.second) {
                    assert(cycle_visited[n_e.first]==cycle_visited[se]);
                    for (auto scb: edge_cycle[edge2id[n_e.first][se]]) {
                        if (!join_visited[scb])cycle_component_cb_set[cycle_visited[se]-1].insert(scb);
                    }
                }
            }
            unordered_set<int> final_candidate=cycle_component_cb_set[0];
            unordered_set<int> other_candidate; other_candidate.clear();
            for (int i=1;i<connect_component;i++)
            {
                for (auto kk:cycle_component_cb_set[0]){
                    other_candidate.insert(kk);
                }
            }
            for (auto kk: other_candidate) {
                if (cycle_component_cb_set[0].find(kk)==cycle_component_cb_set[0].end()) final_candidate.erase(kk);
            }
            if (final_candidate.size()==0) {
                cout << "Terminate" << endl;
                return;
            } else {
                auto nb = -1;
                unordered_set<int> f2=final_candidate;
                for (auto cnb: final_candidate){
                    join_cycles(cnb);
                    if (judge()) {
                        nb=cnb;
                        join_visited[nb]=1;
                        break;
                    }
                    join_cycles(cnb);
                    f2.erase(cnb);
                }
                if (nb==-1) {
                    cout << "Terminate" << endl;
                    return;
                }
                
                unordered_set<int> new_nb; new_nb.clear();
                for (auto eid: current_cir) {
                    for (auto scb: edge_cycle[eid]) {
                        if (!join_visited[scb])new_nb.insert(scb);
                    }
                }
                
                search_full_cycle(1,new_nb);
                join_cycles(nb);
                f2.erase(nb);
                search_full_cycle(0,f2);
                join_visited[nb]=0;
            }
        return;
        }
    }
    if (current_nb.size()==0) return;
    if (!judge()) {
        auto nb = -1;
        unordered_set<int> f2=current_nb;
        for (auto cnb: current_nb){
            join_cycles(cnb);
            if (judge()) {
                nb=cnb;
                join_visited[nb]=1;
                break;
            }
            join_cycles(cnb);
            f2.erase(cnb);
        }
        if (nb==-1) {
            cout << "Terminate" << endl;
            return;
        }
        unordered_set<int> new_nb; new_nb.clear();
        for (auto eid: current_cir) {
            for (auto scb: edge_cycle[eid]) {
                if (!join_visited[scb])new_nb.insert(scb);
            }
        }
        search_full_cycle(1,new_nb);
        join_cycles(nb);
        f2.erase(nb);
        search_full_cycle(0,f2);
        join_visited[nb]=0;
    } else {
        auto nb = *current_nb.begin();
        unordered_set<int> new_nb=current_nb;
        join_cycles(nb);
        join_visited[nb]=1;
        for (auto nb_nb:join_graph[nb]){
            if(join_visited[nb_nb]==0){
                new_nb.insert(nb_nb);
            }
        };
        new_nb.erase(nb);
        search_full_cycle(1,new_nb);
        join_cycles(nb);
        new_nb=current_nb;
        new_nb.erase(nb);
        search_full_cycle(0,new_nb);
        join_visited[nb]=0;
    }
    
}

void search_start_cycle(int index) {
    if (index == start_cb.size()) {
        if (start_cycle.size() % 2 != 0) {
            join_visited.clear();join_visited.resize(cb_num);
            for (auto cid: start_cb) {
                join_visited[cid] = 1;
            }
            current_cir.clear();
            unordered_set<int> current_nb; current_nb.clear();
            for (int cid : start_cycle) {
                join_cycles(cid);
                for (auto nb: join_graph[cid]) {
                    if (join_visited[nb]==0) {
                        current_nb.insert(nb);
                    }
                }
            }
            search_full_cycle(1,current_nb);
        }
        return;
    }

    search_start_cycle(index + 1);

    start_cycle.push_back(start_cb[index]);
    search_start_cycle(index + 1);
    start_cycle.pop_back();
}

string filename, dirname;
int N;

int main(int argc, char **argv) {
    filename = argv[1];
    dirname = argv[2];
    ifstream edge_id_map("results/"+dirname+"/mc_"+filename+"-edge_id");
    int x,y,eid;
    while(edge_id_map >> x >> y >> eid){
        edge2id[x][y] = eid; edge2id[y][x] = eid;
        id2edge[eid] = PII(x,y);
    }
    edge_id_map.close();

    cout << "input s and t: " << endl;
    cin >> s >> t;
    edge_st = edge2id[s][t];
    ifstream cycle_basis_input("results/"+dirname+"/mc_"+filename+"-cycle_basis");
    int p=0;
    cycle_basis_input >> m >> cb_num ;
    cycle_basis.resize(cb_num);
    n = m+1-cb_num;
    sup_graph.resize(n);
    edge_cycle.resize(m);
    string line;
    int num_edge=0;
    getline(cycle_basis_input,line);
    while(getline(cycle_basis_input,line)){
        istringstream iss(line);
        int nn;
        while(iss >> nn){
            if (nn==edge_st) start_cb.push_back(p);
            cycle_basis[p].insert(nn);
            edge_cycle[nn].push_back(p);
            sup_graph[id2edge[nn].first].push_back(p);
            sup_graph[id2edge[nn].second].push_back(p);
                num_edge++;
        }
        p++;
    }
    cycle_basis_input.close();
    cout << "start cycles number: " << start_cb.size() << endl;
    double t1 = clock();
    join_graph.resize(cb_num);
    for (int i=0; i<n; i++){
        if (sup_graph[i].size()==1) continue;
        for (auto x: sup_graph[i]) {
            for (auto y: sup_graph[i]) {
                if (x==y) continue;
                join_graph[x].insert(y);
            }
        }
    }
    sup_graph.clear();
    t2 = clock();
    cout << "build join graph: " << (t2 - t1) / CLOCKS_PER_SEC << endl;
   
    start_cycle.clear();
    search_start_cycle(0);
    cout << circuit_cnt << endl;
    return 0;
}