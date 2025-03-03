#include <bits/stdc++.h>
#include "check.h"
using namespace std;
#define PII pair<int, int>

int m, n;
vector<vector<int>> edge;
vector<int> sp_visited, sp_father, sp_depth, lca_visited, lca_father;
vector<vector<int>> sp_tree_i, query_lca, non_tree_i;
unordered_map<int, unordered_map<int,int>> lcas; 

vector<int> low, dfn, in_stack, s, scc, sz;
int dfncnt = 1, scc_num = 0, tp = 0;
vector<int> cycle1, cycle2;
vector<int> found_cycle;
vector<int> ans; 

unordered_map<int, PII> id2edge;
unordered_map<int, unordered_map<int, int>> edge2id;

string filename;

int dfs_tree_h=0;
vector<unordered_map<int,PII>> change_history;

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

void tarjan_search(int u)
{
    in_stack[u] = 1;
    low[u] = dfncnt;
    dfn[u] = dfncnt;
    dfncnt++;
    s[++tp] = u;
    for (auto v : edge[u])
    {
        if (!dfn[v])
        {
            tarjan_search(v);
            low[u] = min(low[u], low[v]);
        }
        else if (in_stack[v])
        {
            low[u] = min(low[u], dfn[v]);
        }
    }
    if (dfn[u] == low[u])
    {
        ++scc_num;
        while (s[tp] != u)
        {
            scc[s[tp]] = scc_num;
            sz[scc_num]++;
            in_stack[s[tp]] = 0;
            --tp;
        }
        scc[s[tp]] = scc_num;
        sz[scc_num]++;
        in_stack[s[tp]] = 0;
        --tp;
    }
}

void find_scc()
{
    low.resize(n);
    dfn.resize(n);
    in_stack.resize(n);
    s.resize(n+5);
    scc.resize(n);
    sz.resize(n);
    for (int i = 0; i < n; i++)
    {
        if (!dfn[i]){
            tarjan_search(i);
        }
    }
}

struct ST{
    int Max, edge_id;
};
vector<ST> tree;

PII query(int lq,int rq,int l,int r,int x)
{
    if (lq<=l&&rq>=r)
    {
        return PII(tree[x].Max,tree[x].edge_id);
    }
    int mid=l+r>>1,max_gap=0,max_edge_id;
    if (lq<=mid)
    {
        PII tmp=query(lq,rq,l,mid,x<<1);
        if (tmp.first>max_gap) max_gap=tmp.first,max_edge_id=tmp.second;
    }
    if (rq>mid)
    {
        PII tmp=query(lq,rq,mid+1,r,x<<1|1);
        if (tmp.first>max_gap) max_gap=tmp.first,max_edge_id=tmp.second;
    }
    return PII(max_gap,max_edge_id);
}

void update(int q,int gap,int edge_id, int l,int r,int x)
{
    if (l==r)
    {
        tree[x].Max=gap;
        tree[x].edge_id=edge_id;
        return ;
    }
    int mid=l+r>>1;
    if (q<=mid) update(q,gap,edge_id,l,mid,x<<1);
    else update(q,gap,edge_id,mid+1,r,x<<1|1);
    if (tree[x<<1].Max>tree[x<<1|1].Max)
    {
        tree[x].Max=tree[x<<1].Max;
        tree[x].edge_id=tree[x<<1].edge_id;
    }
    else{
        tree[x].Max=tree[x<<1|1].Max;
        tree[x].edge_id=tree[x<<1|1].edge_id;
    }
}

void dfs_search_sptree(int root) {
    for (auto root_nb: edge[root]) {
        if (sp_father[root]==root_nb) continue;
        if (sp_depth[root_nb]<sp_depth[root]) {
            int current_edge = edge2id[root_nb][root];
            int select_edge = -1;
            PII max_gap_pair = query(sp_depth[root_nb]+1,sp_depth[root],1,dfs_tree_h,1);
            if (max_gap_pair.first!=0) select_edge=max_gap_pair.second;
            change_history[root][root_nb]=query(sp_depth[root_nb]+1,sp_depth[root_nb]+1,1,dfs_tree_h,1);
            update(sp_depth[root_nb]+1,sp_depth[root]-sp_depth[root_nb],current_edge,1,dfs_tree_h,1);
        }
    }
    for (auto sp_son: sp_tree_i[root]) {
        dfs_search_sptree(sp_son);
    }
    for (auto pp: change_history[root]){
        update(sp_depth[pp.first]+1,pp.second.first,pp.second.second,1,dfs_tree_h,1);
    }
    return;
}

void dfs(int root,int depth)
{
    sp_visited[root] = 1; sp_depth[root] = depth;
    int ff=false;
    for (auto x: edge[root]) {
        if (sp_visited[x]!=1) {
            ff=true;
            sp_father[x] = root;
            dfs(x,depth+1);
        }
    }
    if (!ff) dfs_tree_h=max(dfs_tree_h,depth+1);
}

void generate_sp_tree(int r) {
    sp_visited.clear(); sp_father.clear(); sp_depth.clear();
    sp_visited.resize(n);
    sp_father.resize(n);
    sp_depth.resize(n);
    sp_father[r]=-1;
    dfs(r,0);
}

bool edge_cmp(int x, int y){
    return sp_depth[x]>sp_depth[y];
}

void cycle_basis() {
    int para=0;
    int set_root = min(para,int(n-1));
    generate_sp_tree(set_root);
    tree.resize(4*dfs_tree_h+100);
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
        sort(edge[x].begin(),edge[x].end(),edge_cmp);
    }

    change_history.clear(); change_history.resize(n);
    dfs_search_sptree(set_root);
}

int main(int argc, char **argv)
{
    // preprocessed graph
    double t1=clock();
    read_graph(argv[1]); 
    filename=argv[1];
    double t2=clock();
    cout << "Reading Graph: " << (t2 - t1) / CLOCKS_PER_SEC << endl;

    int edge_num_id=0;
    for (int u = 0; u < n; u++)
    {
        for (auto v : edge[u])
        {
            if (u>=v) continue;
            edge2id[u][v] = edge_num_id; edge2id[v][u] = edge_num_id;
            id2edge[edge_num_id] = PII(u, v);
            edge_num_id++;
        }
    }
    assert(edge_num_id==m);
    double t5=clock();
    cout << "Edge Serialize: " << (t5 - t2) / CLOCKS_PER_SEC << endl;

    cycle_basis();
    double t6=clock();
    cout << "Cycle Basis: " << (t6 - t5) / CLOCKS_PER_SEC << endl;

    return 0;
}