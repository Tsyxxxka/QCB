#include <bits/stdc++.h>
using namespace std;

int m,n;
unordered_map<int, int> id2name, name2id;
unordered_map<int, unordered_set<int>> graph,graph_final;
vector<vector<int>> edge;
vector<int> visited;

void read_graph(char* dataset){
    char input[100], output[100];
    sprintf(input, "../data/%s", dataset);

    ifstream in(input);
    n = 0;
    m = 0;
    int x, y;
    while (in >> x >> y)
    {
        if (name2id.find(x) == name2id.end())
            name2id[x] = n, id2name[n] = x, n++;
        if (name2id.find(y) == name2id.end())
            name2id[y] = n, id2name[n] = y, n++;
        x = name2id[x];
        y = name2id[y];
        graph[x].insert(y);
        graph[y].insert(x);
    }
    edge.resize(n);
    for (int x = 0; x < n; x++) {
        for (auto y : graph[x]) {
            if (y>=x) continue;
            edge[x].push_back(y), m++;
            edge[y].push_back(x);
        }
    }
    graph.clear(); id2name.clear(); name2id.clear();
    cout << "#vertices: " << n << "\t#edges:" << m << endl;
}

void dfs(int r, int flag) {
    visited[r]=flag;
    for (auto nn: edge[r]) {
        if (!visited[nn]) dfs(nn, flag);
    }
}

int main(int argc, char **argv) {
    read_graph(argv[1]);
    string filename(argv[1]);
    visited.clear(); visited.resize(n);
    int cnt_c=0, mc_s=0, mc=-1;
    for (int i=0; i<n; i++) {
        if (!visited[i]) {
            cnt_c++;
            dfs(i,cnt_c);
            int k=0;
            for (int j=0; j<n; j++) {
                if (visited[j]==cnt_c) k++;
            }
            if (k>mc_s) mc_s=k, mc=cnt_c;
        }
    }
    cout << mc_s << endl;
    int new_n=0;
    for (int i=0; i<n; i++) {
        if (visited[i]!=mc) continue;
        for (auto j: edge[i]) {
            if (visited[j]!=mc) continue;
            if (i>=j) continue;
            if (name2id.find(i) == name2id.end())
                name2id[i] = new_n, id2name[new_n] = i, new_n++;
            if (name2id.find(j) == name2id.end())
                name2id[j] = new_n, id2name[new_n] = j, new_n++;
            graph_final[name2id[i]].insert(name2id[j]);
        }
    }
    assert(new_n==mc_s);
    ofstream outfile("./data/mc_"+filename);
    outfile << mc_s << endl;
    for (int pp=0; pp<mc_s; pp++) {
        for (auto vv: graph_final[pp]) {
            outfile << pp << " " << vv << endl;
        }
    }
    outfile.close();

    
    // ofstream outfile("./data/"+filename);
    // outfile << n << " " << m << endl;
    // for (int i=0; i<n; i++) {
    //     for (auto j: edge[i])
    //         outfile << i << " " << j << endl;
    // }
    // outfile.close();
    return 0;
}