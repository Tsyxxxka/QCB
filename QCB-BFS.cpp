#include <bits/stdc++.h>
#include "check.h"
using namespace std;
#define PII pair<int, int>

int m, n;
vector<vector<int>> edge;
vector<int> sp_visited, sp_father, lca_visited, lca_father;
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
unordered_set<int> bidirected_edge;
string filename;
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

vector<int> find_tree_path(int p, int s) {
    // including p and s, from p to s
    ans.clear(); ans.push_back(s);
    while(s!=p){
        s=sp_father[s];
        if(s!=p)ans.push_back(s);
    }
    ans.push_back(p);
    reverse(ans.begin(),ans.end());
    return ans;
} 

void search_cycle(int u, int v)
{
    int lca_uv = lcas[u][v];
    if (lca_uv==u || lca_uv==v) {
        int p1=lca_uv, p2=lca_uv==u?v:u;
        vector<int> uv_path = find_tree_path(p1,p2);
        assert(uv_path.size()>=3);
        cycle1.clear(); 
        cycle1.push_back(edge2id[u][v]);
        int path_len=uv_path.size()-1;
        int cs=path_len-1, flag=1;
        while(flag && cs>1){
            for (int n1=0; n1+cs<=path_len; n1++){
                int pp = uv_path[n1], ss = uv_path[n1+cs];
                if (!edge[pp].empty() && find(edge[pp].begin(),edge[pp].end(),ss)!=edge[pp].end()){
                    for (int k=0; k<n1; k++) cycle1.push_back(edge2id[uv_path[k]][uv_path[k+1]]);
                    cycle1.push_back(edge2id[pp][ss]);
                    for (int k=n1+cs; k<path_len; k++) cycle1.push_back(edge2id[uv_path[k]][uv_path[k+1]]);
                    flag=0;
                    break;
                }   
            }
            cs--;
        }
        if (flag==1) {
            for (int k=0; k<path_len; k++) cycle1.push_back(edge2id[uv_path[k]][uv_path[k+1]]);
        }
    } else {
        cycle1.clear(); cycle2.clear();
        cycle1.push_back(edge2id[u][v]);
        int flag=1, p1=u, p2=v;
        while(flag) {
            cycle2 = cycle1;
            p2=v;
            while(p2!=lca_uv) {
                if ((p1!=u || p2!=v ) && !edge[p2].empty() && find(edge[p2].begin(),edge[p2].end(),p1)!=edge[p2].end() ) {
                    cycle2.push_back(edge2id[p1][p2]);
                    flag=0;
                    break;
                }
                if (sp_father[p2]!=lca_uv) {
                    cycle2.push_back(edge2id[p2][sp_father[p2]]);
                    p2=sp_father[p2];
                } else break;
            }
            if (flag==0) break;
            if (sp_father[p1]!=lca_uv) {
                cycle1.push_back(edge2id[p1][sp_father[p1]]);
                p1=sp_father[p1];
            } else break;
        }
        if (flag==0);
        else {
            cycle2.push_back(edge2id[p1][lca_uv]);
            cycle2.push_back(edge2id[p2][lca_uv]);
        }
    }
}


void dfs(int root)
{
    sp_visited[root] = 1;
    for (auto x: edge[root]) {
        if (sp_visited[x]!=1) {
            sp_father[x] = root;
            dfs(x);
        }
    }
}

void bfs(int root) {
    sp_visited[root] = 1;
    queue<int> Q;
    Q.push(root);
    while(!Q.empty()) {
        int f = Q.front();
        Q.pop();
        for (auto ff: edge[f]){
            if (sp_visited[ff]!=1) sp_visited[ff]=1,sp_father[ff]=f, Q.push(ff);
        }
    }
}

void generate_sp_tree(int r) {
    sp_visited.clear(); sp_father.clear();
    sp_visited.resize(n);
    sp_father.resize(n);
    sp_father[r]=-1;
    bfs(r);
}

int find(int q) {
    if (lca_father[q]==q) return q;
    else return lca_father[q]=find(lca_father[q]);
}

void dfs_lca(int root) {
    lca_father[root] = root; lca_visited[root] = 1;
    for (auto x: sp_tree_i[root]) {
        if (lca_visited[x]!=1) {
            dfs_lca(x);
            lca_father[x]=root;
        }
    }
    for (auto q: query_lca[root]) {
        if (lca_visited[q] && lcas[q][root]==-1) {
            int lca_ans = find(q);
            lcas[q][root] = lca_ans; lcas[root][q] = lca_ans;
        }
    }
}

void find_LCA(int r) {
    lca_visited.clear(); lca_father.clear();
    lca_visited.resize(n,0);
    lca_father.resize(n,0);
    dfs_lca(r);
}

void cycle_basis() {
    int para=0;
    int set_root = min(para,int(n-1));
    generate_sp_tree(0);
    sp_tree_i.clear(); query_lca.clear(); non_tree_i.clear();
    sp_tree_i.resize(n); query_lca.resize(n); non_tree_i.resize(n);
    lcas.clear();
    for (int x=0; x<n; x++) {
        for (auto y: edge[x]) {
            if (y<=x) continue;
            if (sp_father[x]==y) {
                sp_tree_i[y].push_back(x);
            } else if (sp_father[y]==x) {
                sp_tree_i[x].push_back(y);
            } else {
                lcas[x][y] = -1; lcas[y][x] = -1;
                query_lca[x].push_back(y); query_lca[y].push_back(x);
                non_tree_i[x].push_back(y);
            }
        }
    }
    find_LCA(set_root);

    edge_cnt.resize(m);
    outfile.open("results/QCB_BFS/"+filename+"-cycle_basis");
    outfile << m << " " << m-n+1 << endl;
    for (int x=0; x<n; x++) {
        for (auto y: non_tree_i[x]) {
            search_cycle(x, y);
            for (auto x: cycle2) {
                edge_cnt[x]++;
                outfile << x << " ";
            }
            outfile << endl;
        }
    }
    outfile.close();
}

int main(int argc, char **argv)
{
    // preprocessed graph
    double t1=clock();
    read_graph(argv[1]); 
    filename = argv[1];
    double t2=clock();
    cout << "Reading Graph: " << (t2 - t1) / CLOCKS_PER_SEC << endl;

    ofstream edge_id_map("results/QCB_BFS/"+filename+"-edge_id");
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

    // int max_edge_cnt=0, total_cnt=0;
    // for (int i=0; i<m; i++) {
    //     total_cnt+=edge_cnt[i];
    //     max_edge_cnt=max(max_edge_cnt,edge_cnt[i]);
    // }
    // cout << "MAX EDGE CNT: " << max_edge_cnt << endl;
    // cout << "AVG EDGE CNT: " << total_cnt*1.0/m << endl;

    if (argv[2]) {
        vector<int> zeroRow = check_rank(m,m-n+1,4,filename);
        cout << "Zero row number: " << zeroRow.size() << endl;
    }

    return 0;
}