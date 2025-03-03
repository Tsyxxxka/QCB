#include <bits/stdc++.h>
#include "check.h"
using namespace std;
#define PII pair<int, int>

int m, n;
vector<vector<int>> edge;
vector<int> sp_visited, sp_father, sp_depth, lca_visited, lca_father;
vector<vector<int>> sp_tree_i, query_lca, non_tree_i;
unordered_map<int, unordered_map<int,int>> lcas; 
// int reduced_edge = 0;
vector<int> edge_cnt;
vector<int> son_cnt;

// scc[i]=k: node i in k-th scc (counting from 1)
// sz[k]=s: k-th scc has s nodes
vector<int> low, dfn, in_stack, s, scc, sz;
int dfncnt = 1, scc_num = 0, tp = 0;
vector<int> cycle1, cycle2;
vector<int> found_cycle;
vector<int> ans; 

unordered_map<int, PII> id2edge;
unordered_map<int, unordered_map<int, int>> edge2id;

unordered_set<int> bidirected_edge;
string filename;

class cmp {
public:
    bool operator ()(const PII &a,const PII &b) {
        return a.second>b.second;
    }
};

// set<PII,cmp> checked_non_tree_edge;
int dfs_tree_h=0;
vector<PII> up_nodes_current_max_gap;
vector<unordered_map<int,PII>> change_history;


// struct SCC_Graph{
//     int node_num;
//     int edge_num;
//     // nodes[i]=j: i-th node in current subgraph is the j-th node in original graph
//     vector<int> nodes;
//     // using current node id
//     vector<unordered_set<int>> edges;
//     // g2l_nodes[i] = j : i-th node in original graph is the j-th node in current graph
//     unordered_map<int,int> g2l_nodes;
//     SCC_Graph(): node_num(), edge_num(), nodes(), edges(){};
//     SCC_Graph(int n, int m, vector<int> ns, vector<unordered_set<int>> e): node_num(n), edge_num(m), nodes(ns), edges(e){};
// };

// vector<SCC_Graph> subG;
// vector<unordered_set<int>> tmpE;
// vector<int> node_l2g;
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
            // father[v] = u;
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
    // father.resize(n);
    for (int i = 0; i < n; i++)
    {
        if (!dfn[i]){
            // father[i] = -1
            tarjan_search(i);
        }
    }
}

// vector<int> find_tree_path(int p, int s){
//     // excluding p and s, from p to s
//     vector<int> ans; ans.clear();
//     while (s!=p) {
//         s=sp_father[s];
//         if (s!=p) ans.push_back(s);
//     }
//     reverse(ans.begin(), ans.end());
//     return ans;
// }
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

struct ST{
    int Max, edge_id;
};
vector<ST> tree;

PII query(int lq,int rq,int l,int r,int x)//query[lq,rq] query(lq,rq,1,h,1)
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
    // cout << root << ": ";
    // for (int i=0; i<sp_depth[root]; i++) {
    //     cout << up_nodes_current_max_gap[i].first << up_nodes_current_max_gap[i].second << endl;
    // }
    for (auto root_nb: edge[root]) {
        // cout << root_nb << endl;
        if (sp_father[root]==root_nb) continue;
        if (sp_depth[root_nb]<sp_depth[root]) {
            // cout << root_nb << endl;
            int current_edge = edge2id[root_nb][root];
            // search for max_gap
            // int select_edge = -1, max_gap = 0;
            // for (int ii=sp_depth[root_nb]; ii<sp_depth[root]; ii++){
            //     if (up_nodes_current_max_gap[ii].first>max_gap) {
            //         max_gap = up_nodes_current_max_gap[ii].first;
            //         select_edge = up_nodes_current_max_gap[ii].second;
            //     }
            // }
            int select_edge = -1;
            PII max_gap_pair = query(sp_depth[root_nb]+1,sp_depth[root],1,dfs_tree_h,1);
            if (max_gap_pair.first!=0) select_edge=max_gap_pair.second;
            // if (root_nb==0) {
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
            // }
            // // if (sp_depth[root]-sp_depth[root_nb]>up_nodes_current_max_gap[sp_depth[root_nb]].first){
            // change_history[root][root_nb]=up_nodes_current_max_gap[sp_depth[root_nb]];
            // up_nodes_current_max_gap[sp_depth[root_nb]]=PII(sp_depth[root]-sp_depth[root_nb],current_edge);
            // // }
            // change_history[i][j]=PPI(max_gap, max_edge_id): dfs_search to i makes i's ancestor j's max_gap changes from (max_gap, max_edge_id) to (new_gap, edge_ij)

            change_history[root][root_nb]=query(sp_depth[root_nb]+1,sp_depth[root_nb]+1,1,dfs_tree_h,1);
            update(sp_depth[root_nb]+1,sp_depth[root]-sp_depth[root_nb],current_edge,1,dfs_tree_h,1);
        }
    }
    for (auto sp_son: sp_tree_i[root]) {
        // up_nodes_current_max_gap.push_back(PII(-1,-1));
        dfs_search_sptree(sp_son);
        // up_nodes_current_max_gap.pop_back();
    }
    for (auto pp: change_history[root]){
        // up_nodes_current_max_gap[sp_depth[pp.first]]=pp.second;
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

// get sp_father[]
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
    // cout << "DFS_TREE_H: " << dfs_tree_h << endl;
    tree.resize(4*dfs_tree_h+100);
    // for (int i=0; i<n; i++) cout << i << " " << sp_depth[i]<< endl;
    sp_tree_i.clear();
    sp_tree_i.resize(n);
    for (int x=0; x<n; x++) {
        for (auto y: edge[x]) {
            // for undirected graph subG, only record once for each undirected edge
            if (y<=x) continue;
            if (sp_father[x]==y) {
                sp_tree_i[y].push_back(x);
            } else if (sp_father[y]==x) {
                sp_tree_i[x].push_back(y);
            }
        }
        sort(edge[x].begin(),edge[x].end(),edge_cmp);
    }

    // for (int x=0; x<n; x++) {
    //     cout << x << ": ";
    //     for (auto y: edge[x]) {
    //         cout << y << " ";
    //     }
    //     cout << endl;
    // }


    edge_cnt.resize(m);
    outfile << m << " " << m-n+1 << endl;
    // checked_non_tree_edge.clear();
    // up_nodes_current_max_gap.clear();
    change_history.clear(); change_history.resize(n);
    dfs_search_sptree(set_root);
    outfile.close();
    // for (int x=0; x<n; x++) {
    //     for (auto y: non_tree_i[x]) {
    //         search_cycle(x, y);
    //         for (auto x: cycle1) {
    //             edge_cnt[x]++;
    //             outfile << x << " ";
    //         }
    //         outfile << endl;
    //     }
    // }
}

int main(int argc, char **argv)
{
    // underlying undirected graph should be connected
    double t1=clock();
    read_graph(argv[1]); 
    filename=argv[1];
    cout << filename << endl;
    outfile.open("results/ST/"+filename+"-cycle_basis");
    double t2=clock();
    cout << "Reading Graph: " << (t2 - t1) / CLOCKS_PER_SEC << endl;

    ofstream edge_id_map("results/ST/"+filename+"-edge_id");
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
        vector<int> zeroRow = check_rank(m,m-n+1,1,filename);
        cout << "Zero row number: " << zeroRow.size() << endl;
    }

    return 0;
}