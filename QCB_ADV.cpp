#include <bits/stdc++.h>
#include "check.h"
using namespace std;
#define PII pair<int, int>

int m, n;
vector<vector<int>> edge;
vector<int> sp_visited, sp_father, sp_depth, lca_visited, lca_father;
vector<vector<int>> sp_tree_i, query_lca, non_tree_i;
unordered_map<int, unordered_map<int,int>> lcas;
vector<int> edge_cnt;
vector<int> son_cnt;

vector<int> low, dfn, in_stack, s, scc, sz;
int dfncnt = 1, scc_num = 0, tp = 0;
vector<int> cycle1, cycle2;
vector<int> found_cycle;
vector<int> ans; 

unordered_map<int, PII> id2edge;
unordered_map<int, unordered_map<int, int>> edge2id;

string filename;

vector<PII> up_nodes_current_max_gap;
vector<unordered_map<int,PII>> change_history;

ofstream outfile;

void read_graph(char *dataset)
{
    char input[100], output[100];
    sprintf(input, "./data/%s", dataset);

    ifstream in(input);
    in >> n;
    edge.resize(n);
    int x, y;
    while (in >> x >> y) {
        m++;
        edge[x].push_back(y);
        edge[y].push_back(x);
    }
    cout << "#vertices: " << n << "\t#edges:" << m << endl;
}

void find_tree_path(int p, int s) {
    // from p to s, all edge id
    ans.clear(); 
    if (s==p) return;
    while (sp_father[s]!=p) {
        ans.push_back(edge2id[s][sp_father[s]]);
        s=sp_father[s];
    }
    ans.push_back(edge2id[s][p]);
    reverse(ans.begin(),ans.end());
} 

void dfs_search_sptree(int root) {
    for (auto root_nb: edge[root]) {
        if (sp_father[root]==root_nb) continue;
        if (sp_depth[root_nb]<sp_depth[root]) {
            int current_edge = edge2id[root_nb][root];
            int select_edge = -1, max_gap = 0;
            for (int ii=sp_depth[root_nb]; ii<sp_depth[root]; ii++){
                if (up_nodes_current_max_gap[ii].first>max_gap) {
                    max_gap = up_nodes_current_max_gap[ii].first;
                    select_edge = up_nodes_current_max_gap[ii].second;
                }
            }
            outfile << current_edge << " ";
            if (select_edge==-1) {
                find_tree_path(root_nb, root); 
                for (auto xx: ans) outfile << xx << " ";
                outfile << endl;
            } else {
                int acs_node = id2edge[select_edge].first, dec_node = id2edge[select_edge].second;
                if (sp_depth[acs_node]>sp_depth[dec_node]) swap(acs_node, dec_node);
                find_tree_path(root_nb,acs_node);
                for (auto xx: ans) outfile << xx << " ";
                outfile << select_edge << " ";
                find_tree_path(dec_node,root);
                for (auto xx: ans) outfile << xx << " ";
                outfile << endl;
            }
            if (sp_depth[root]-sp_depth[root_nb]>up_nodes_current_max_gap[sp_depth[root_nb]].first){
                change_history[root][root_nb]=up_nodes_current_max_gap[sp_depth[root_nb]];
                up_nodes_current_max_gap[sp_depth[root_nb]]=PII(sp_depth[root]-sp_depth[root_nb],current_edge);
            }
        }
    }
    for (auto sp_son: sp_tree_i[root]) {
        up_nodes_current_max_gap.push_back(PII(-1,-1));
        dfs_search_sptree(sp_son);
        up_nodes_current_max_gap.pop_back();
    }
    for (auto pp: change_history[root]){
        up_nodes_current_max_gap[sp_depth[pp.first]]=pp.second;
    }
    return;
}


void dfs(int root,int depth)
{
    sp_visited[root] = 1; sp_depth[root] = depth;
    for (auto x: edge[root]) {
        if (sp_visited[x]!=1) {
            sp_father[x] = root;
            dfs(x,depth+1);
        }
    }
}

void generate_sp_tree(int r) {
    sp_visited.clear(); sp_father.clear(); sp_depth.clear();
    sp_visited.resize(n);
    sp_father.resize(n);
    sp_depth.resize(n);
    sp_father[r]=-1;
    dfs(r,0);
}

void cycle_basis() {
    int para=0;
    int set_root = min(para,int(n-1));
    generate_sp_tree(set_root);
    sp_tree_i.clear();
    sp_tree_i.resize(n);
    for (int x=0; x<n; x++) {
        for (auto y: edge[x]) {
            if (y<=x) continue;
            if (sp_father[x]==y) {
                sp_tree_i[y].push_back(x);
            } else if (sp_father[y]==x) {
                sp_tree_i[x].push_back(y);
            }
        }
    }

    edge_cnt.resize(m);
    
    outfile << m << " " << m-n+1 << endl;
    up_nodes_current_max_gap.clear();
    change_history.clear(); change_history.resize(n);
    dfs_search_sptree(set_root);
    outfile.close();
}

int main(int argc, char **argv)
{
    // preprocessed graph
    double t1=clock();
    read_graph(argv[1]); 
    filename=argv[1];
    outfile.open("results/ADV/"+filename+"-cycle_basis");
    double t2=clock();
    cout << "Reading Graph: " << (t2 - t1) / CLOCKS_PER_SEC << endl;

    ofstream edge_id_map("results/ADV/"+filename+"-edge_id");
    int edge_num_id=0;
    for (int u = 0; u < n; u++)
    {
        for (auto v : edge[u])
        {
            if (u>=v) continue;
            edge2id[u][v] = edge_num_id; edge2id[v][u] = edge_num_id;
            id2edge[edge_num_id] = PII(u, v);
            edge_id_map << u << " " << v << " " << edge_num_id << endl;
            edge_num_id++;
        }
    }
    assert(edge_num_id==m);
    edge_id_map.close();
    double t5=clock();
    cout << "Edge Serialize: " << (t5 - t2) / CLOCKS_PER_SEC << endl;

    cycle_basis();
    double t6=clock();
    cout << "Cycle Basis: " << (t6 - t5) / CLOCKS_PER_SEC << endl;

    if (argv[2]) {
        vector<int> zeroRow = check_rank(m,m-n+1,2,filename);
        cout << "Zero row number: " << zeroRow.size() << endl;
    }

    return 0;
}