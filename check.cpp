#include "check.h"
#include <bits/stdc++.h>

using namespace std;

const double EPS = 1e-9;

int gaussianRank(vector<vector<double>> &matrix, vector<int> &zeroRows) {
    int m = matrix.size();
    int n = matrix[0].size();
    int rank = 0;

    vector<int> rowIndices(m);
    for (int i = 0; i < m; ++i) {
        rowIndices[i] = i;
    }

    for (int col = 0, row = 0; col < n && row < m; ++col) {
        int pivot = row;
        for (int i = row + 1; i < m; ++i) {
            if (fabs(matrix[i][col]) > fabs(matrix[pivot][col])) {
                pivot = i;
            }
        }

        if (fabs(matrix[pivot][col]) < EPS)
            continue;

        swap(matrix[row], matrix[pivot]);
        swap(rowIndices[row], rowIndices[pivot]);

        double pivotValue = matrix[row][col];
        for (int j = col; j < n; ++j) {
            matrix[row][j] /= pivotValue;
        }

        for (int i = 0; i < m; ++i) {
            if (i != row && fabs(matrix[i][col]) > EPS) {
                double factor = matrix[i][col];
                for (int j = col; j < n; ++j) {
                    matrix[i][j] -= factor * matrix[row][j];
                }
            }
        }
        ++rank;
        ++row;
    }

    for (int i = 0; i < m; ++i) {
        bool isZeroRow = true;
        for (int j = 0; j < n; ++j) {
            if (fabs(matrix[i][j]) > EPS) {
                isZeroRow = false;
                break;
            }
        }
        if (isZeroRow) {
            zeroRows.push_back(rowIndices[i]);
        }
    }

    return rank;
}

vector<int> check_rank(int m, int n, int t, string filename) {
    string dirname;
    if (t==1) {
        dirname = "ST/";
    } else if (t==2) {
        dirname = "ADV/";
    } else if (t==3) {
        dirname = "BAC/";
    } else {
        dirname = "QCB_BFS/";
    }
    ifstream in("./results/"+dirname+filename+"-cycle_basis");
    cout << m << " " << n << endl;
    vector<vector<double>> matrix(m, vector<double>(n));
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < m; ++i) {
            matrix[i][j]=0;
        }
    }
    string line;
    getline(in,line);
    int col=0;
    while(getline(in,line)){
        istringstream iss(line);
        int nn;
        while(iss >> nn){
            matrix[nn][col]=1;
        }
        col++;
    }
    
    vector<int> zeroRows;
    int rank = gaussianRank(matrix, zeroRows);

    cout << "Rank: " << rank << endl;
    if (rank==n) cout << "GOOD" << endl;
    else cout << "ERROR" << endl;

    return zeroRows;
}

