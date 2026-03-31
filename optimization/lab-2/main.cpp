#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>
#include <set>
#include <cmath>
#include <utility>

using namespace std;

const double EPS = 1e-10;

// ------------------------------------------------------------
// Construction of initial basic plan using the north-west corner method
// with degeneracy handling (adding zero basic cells)
// ------------------------------------------------------------
void northWestAngle(const vector<double>& a, const vector<double>& b,
                    vector<vector<double>>& x, vector<pair<int,int>>& basis) {
    int m = a.size();
    int n = b.size();
    x.assign(m, vector<double>(n, 0.0));
    basis.clear();

    vector<double> supply = a;
    vector<double> demand = b;

    int i = 0, j = 0;
    while (i < m && j < n) {

        if (supply[i] < EPS) {
            i++;
            continue;
        }
        if (demand[j] < EPS) {
            j++;
            continue;
        }

        double amount = min(supply[i], demand[j]);
        x[i][j] = amount;
        basis.push_back({i, j});

        supply[i] -= amount;
        demand[j] -= amount;

        // degeneracy
        if (supply[i] < EPS && demand[j] < EPS) {
            
            if (i + 1 < m) {
         
                basis.push_back({i + 1, j});    
                i++;                           
            } else if (j + 1 < n) {
          
                basis.push_back({i, j + 1});    
                j++;                            
            }
           
        }
        else if (supply[i] < EPS) {
            i++;
        }
        else if (demand[j] < EPS) {
            j++;
        }
    }


    // if ((int)basis.size() != m + n - 1) {
    //     cerr << "Warning: basis size = " << basis.size()
    //          << ", required " << m + n - 1 << ". Possible problems.\n";
        
    // }
}

// ------------------------------------------------------------
// Calculation of potentials u and v from basic cells
// ------------------------------------------------------------
bool computePotentials(const vector<vector<double>>& c,
                       const vector<pair<int,int>>& basis,
                       vector<double>& u, vector<double>& v) {
    int m = c.size();
    int n = c[0].size();

    u.assign(m, numeric_limits<double>::quiet_NaN());
    v.assign(n, numeric_limits<double>::quiet_NaN());

    u[0] = 0.0;  

    bool changed;
    do {
        changed = false;
        for (const auto& cell : basis) {
            int i = cell.first;
            int j = cell.second;
            if (!isnan(u[i]) && isnan(v[j])) {
                v[j] = c[i][j] - u[i];
                changed = true;
            }
            else if (!isnan(v[j]) && isnan(u[i])) {
                u[i] = c[i][j] - v[j];
                changed = true;
            }
        }
    } while (changed);

    
    for (int i = 0; i < m; ++i)
        if (isnan(u[i])) return false;
    for (int j = 0; j < n; ++j)
        if (isnan(v[j])) return false;

    return true;
}

// ------------------------------------------------------------
// Search for a cell with negative evaluation (entering the basis)
// ------------------------------------------------------------
pair<int,int> findEnteringCell(const vector<vector<double>>& c,
                               const vector<double>& u,
                               const vector<double>& v,
                               const vector<pair<int,int>>& basis,
                               double& minDelta) {
    int m = c.size();
    int n = c[0].size();
    set<pair<int,int>> basisSet(basis.begin(), basis.end());

    minDelta = 0.0;
    pair<int,int> entering = {-1, -1};

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            if (basisSet.count({i, j})) continue;
            double delta = c[i][j] - (u[i] + v[j]);
            if (delta < -EPS) {
                if (entering.first == -1 || delta < minDelta) {
                    minDelta = delta;
                    entering = {i, j};
                }
            }
        }
    }
    return entering;
}

// ------------------------------------------------------------
// Construction of the recalculation cycle for the entering cell (recursive depth-first search)
// ------------------------------------------------------------
bool buildCycleDFS(pair<int,int> current,
                   pair<int,int> start,
                   int lastDir,            
                   vector<pair<int,int>>& path,
                   const set<pair<int,int>>& available) {
    if (current == start && !path.empty()) {
        return true;  
    }

    if (lastDir != 0) {
        for (const auto& cell : available) {
            if (cell.first == current.first && cell.second != current.second) {
                if (cell == start && !path.empty()) {
                    path.push_back(cell);
                    return true;
                }
                if (find(path.begin(), path.end(), cell) == path.end()) {
                    path.push_back(cell);
                    if (buildCycleDFS(cell, start, 0, path, available))
                        return true;
                    path.pop_back();
                }
            }
        }
    }

   
    if (lastDir != 1) {
        for (const auto& cell : available) {
            if (cell.second == current.second && cell.first != current.first) {
                if (cell == start && !path.empty()) {
                    path.push_back(cell);
                    return true;
                }
                if (find(path.begin(), path.end(), cell) == path.end()) {
                    path.push_back(cell);
                    if (buildCycleDFS(cell, start, 1, path, available))
                        return true;
                    path.pop_back();
                }
            }
        }
    }

    return false;
}

// Wrapper for cycle construction
vector<pair<int,int>> buildCycle(pair<int,int> entering,
                                 const vector<pair<int,int>>& basis) {
    set<pair<int,int>> available(basis.begin(), basis.end());
    available.insert(entering);

    vector<pair<int,int>> path;
    if (buildCycleDFS(entering, entering, -1, path, available)) {
        vector<pair<int,int>> cycle;
        cycle.push_back(entering);
        cycle.insert(cycle.end(), path.begin(), path.end());
        // Remove possible duplicate entering at the end
        if (!cycle.empty() && cycle.back() == entering && cycle.size() > 1) {
            cycle.pop_back();
        }
        return cycle;
    } else {
        cerr << "Error: failed to build recounting cycle!\n";
        return {};
    }
}

// ------------------------------------------------------------
// Determination of the redistribution amount theta
// ------------------------------------------------------------
double findTheta(const vector<pair<int,int>>& cycle,
                 const vector<vector<double>>& x) {
    double theta = numeric_limits<double>::max();
    // Cells with odd index (starting from 1) have sign "-"
    for (size_t k = 1; k < cycle.size(); ++k) {
        if (k % 2 == 1) {
            int i = cycle[k].first;
            int j = cycle[k].second;
            if (x[i][j] < theta)
                theta = x[i][j];
        }
    }
    return theta;
}

// ------------------------------------------------------------
// Redistribution of supplies along the cycle and basis update
// ------------------------------------------------------------
void redistribute(vector<vector<double>>& x,
                  const vector<pair<int,int>>& cycle,
                  double theta,
                  vector<pair<int,int>>& basis,
                  pair<int,int> entering) {
    // Change values in cycle cells
    for (size_t k = 0; k < cycle.size(); ++k) {
        int i = cycle[k].first;
        int j = cycle[k].second;
        if (k % 2 == 0) {   // "+"
            x[i][j] += theta;
        } else {            // "-"
            x[i][j] -= theta;
        }
    }

    // Determine the leaving cell (one of the cells with sign "-", where x became zero)
    pair<int,int> leaving = {-1, -1};
    for (size_t k = 1; k < cycle.size(); k += 2) {
        int i = cycle[k].first;
        int j = cycle[k].second;
        if (x[i][j] < EPS) {
            leaving = {i, j};
            break;
        }
    }

    if (leaving.first == -1) {
        cerr << "Error: unreachable cell found (theta = " << theta << ")\n";
        return;
    }

    auto it = find(basis.begin(), basis.end(), leaving);
    if (it != basis.end())
        basis.erase(it);
    else
        cerr << "Warning: reachable cell is not in basis\n";

    basis.push_back(entering);
}

// ------------------------------------------------------------
// Calculation of total plan cost
// ------------------------------------------------------------
double totalCost(const vector<vector<double>>& x,
                 const vector<vector<double>>& c) {
    double cost = 0.0;
    int m = x.size();
    int n = x[0].size();
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            cost += x[i][j] * c[i][j];
    return cost;
}

// ------------------------------------------------------------
// Main algorithm of the potential method
// ------------------------------------------------------------
void transportProblemSolver(const vector<double>& a,
                            const vector<double>& b,
                            const vector<vector<double>>& c) {
    int m = a.size();
    int n = b.size();

    // Check balance
    double sumA = 0, sumB = 0;
    for (double val : a) sumA += val;
    for (double val : b) sumB += val;
    if (fabs(sumA - sumB) > EPS) {
        cerr << "Error: task is not balanced (Total supplies = " << sumA
             << ", total demands = " << sumB << ")\n";
        return;
    }

    // Step 1: initial plan using north-west corner method
    vector<vector<double>> x;
    vector<pair<int,int>> basis;
    northWestAngle(a, b, x, basis);

    cout << "Initial plan (NWC method):\n";
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << x[i][j] << "\t";
        }
        cout << "\n";
    }
    cout << "Initial cost: " << totalCost(x, c) << "\n\n";
     vector<double> u, v;
    // Iterative process
    int iter = 0;
    while (true) {
        // Step 2: calculation of potentials
       
        if (!computePotentials(c, basis, u, v)) {
            cerr << "Error: failed to compute potential (basis is not conjugated).\n";
            break;
        }

        // Step 3: search for a cell with negative evaluation
        double minDelta;
        pair<int,int> entering = findEnteringCell(c, u, v, basis, minDelta);

        if (entering.first == -1) {
            cout << "Optimal plan achieved!\n";
            break;
        }

        cout << "Iteration " << ++iter << ": entering cell ("
             << entering.first << ", " << entering.second << ") with evaluation " << minDelta << "\n";

        // Step 4: cycle construction
        vector<pair<int,int>> cycle = buildCycle(entering, basis);
        if (cycle.empty()) break;

        // Step 5: determination of theta
        double theta = findTheta(cycle, x);
        cout << "Theta = " << theta << "\n";

        // Step 6: redistribution
        redistribute(x, cycle, theta, basis, entering);

        // Output current plan
        cout << "New plan:\n";
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                cout << x[i][j] << "\t";
            }
            cout << "\n";
        }
        cout << "Current cost: " << totalCost(x, c) << "\n\n";
    }

    // Final result
    cout << "=========================\n";
    cout << "Optimal transportation plan:\n";
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << x[i][j] << "\t";
        }
        cout << "\n";

    }
    
    cout << "Total cost: " << totalCost(x, c) << "\n";

    cout << "d " << "\n";
     for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            double delta = c[i][j] - (u[i] + v[j]);
            cout << delta << "\t";
        }
         cout << "\n";
    }
    cout << "\n";

    for (int i = 0; i < m; ++i) {
         cout << u[i] << "\t";
    }
    cout << "\n";
    for (int j = 0; j < n; ++j) {
        cout << v[j] << "\t";
    }
    cout << "\n";
}


// ------------------------------------------------------------
// Function for unbalanced problem with surplus from suppliers
// ------------------------------------------------------------
void transportRatingProblemSolver(const vector<double>& a,
                                   const vector<double>& b,
                                   const vector<vector<double>>& c,
                                   const vector<double>& r) {
    int m = a.size();               // number of suppliers
    int n = b.size();               // number of consumers

   
    if ((int)r.size() != m) {
        cerr << "Error: rating vector size does not match the number of suppliers.\n";
        return;
    }

    double totalSupply = 0.0, totalDemand = 0.0;
    for (double val : a) totalSupply += val;
    for (double val : b) totalDemand += val;

   
    if (totalSupply < totalDemand - 1e-10) {
        cerr << "Error: total supply is less than total demand. The problem has no solution.\n";
        return;
    }

    if (fabs(totalSupply - totalDemand) < 1e-10) {
        cout << "The problem is already balanced. Solving with standard method.\n";
        transportProblemSolver(a, b, c);
        return;
    }

    // Surplus to be removed
    double surplus = totalSupply - totalDemand;
    cout << "Cargo surplus: " << surplus << "\n";

   
    vector<int> indices(m);
    for (int i = 0; i < m; ++i) indices[i] = i;
    sort(indices.begin(), indices.end(),
         [&r](int i, int j) { return r[i] < r[j]; });

   
    vector<double> a_balanced = a;

    
    for (int idx : indices) {
        if (surplus <= 0) break;          // surplus fully covered
        double take = min(a_balanced[idx], surplus);
        a_balanced[idx] -= take;
        surplus -= take;
        cout << "From supplier " << idx+1 << " (rating " << r[idx] << ") taken " << take
             << " cargo. Remaining: " << a_balanced[idx] << "\n";
    }

    
    if (surplus > 1e-10) {
        cerr << "Error: failed to completely eliminate surplus. Remaining: " << surplus << "\n";
        return;
    }

    cout << "\nBalanced problem:\n";
    cout << "Supplies after adjustment: ";
    for (double val : a_balanced) cout << val << " ";
    cout << "\n\n";

    // Solve the obtained balanced problem with the standard method
    transportProblemSolver(a_balanced, b, c);
}


int main() {
  
   vector<double> a = {6, 11, 5, 5};   // supplies
    vector<double> b = {6, 3, 4, 8, 6}; // demands

    vector<vector<double>> c = {
    {3, 4, 4, 3, 1},
    {5, 6, 5, 3, 1},
    {4, 3, 2, 6, 1},
    {5, 8, 6, 4, 5}
};

    vector<double> b_1 = {5, 2, 3, 7, 5}; // demands
    vector<double> r = {0,3,4,1}; // ratings
    //transportProblemSolver(a, b, c);

    transportRatingProblemSolver(a, b_1, c, r);

    return 0;
}