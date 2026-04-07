#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <string>
#include <limits>

using namespace std;

// ============================================================
// Simplex Method (two-phase, minimisation)
// ============================================================
struct Simplex {
    int m, n;
    vector<int> B, N;
    vector<vector<double>> D;

    static const double EPS_S;
    static const double INF_S;

    static string varName(int idx, int n) {
        if (idx == -1) return "art";
        if (idx < n)   return "x" + to_string(idx + 1);
        return "s" + to_string(idx - n + 1);
    }

    Simplex(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& c) {
        m = (int)b.size();
        n = (int)c.size();
        B.resize(m);
        N.resize(n + 1);
        D.assign(m + 2, vector<double>(n + 2, 0.0));

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                D[i][j] = A[i][j];

        for (int i = 0; i < m; i++) {
            B[i] = n + i;
            D[i][n]     = -1.0;     // artificial variable column
            D[i][n + 1] = b[i];     // RHS
        }
        for (int j = 0; j < n; j++) {
            N[j] = j;
            D[m][j] = -c[j];        // phase-2 objective row (minimise)
        }
        N[n] = -1;
        D[m + 1][n] = 1.0;          // phase-1 objective
    }

    void logBasis(int phase, int iter) const {
        int obj_row = (phase == 1) ? m + 1 : m;
        cout << "  [Phase " << phase << " | iter " << iter << "]  Basis: {";
        for (int i = 0; i < m; i++) {
            cout << varName(B[i], n);
            if (i + 1 < m) cout << ", ";
        }
        cout << "}  obj=" << fixed << setprecision(4) << D[obj_row][n + 1] << "\n";
        cout << "    BFS: ";
        for (int i = 0; i < m; i++)
            cout << varName(B[i], n) << "=" << fixed << setprecision(2) << D[i][n + 1] << "  ";
        cout << "\n";
    }

    void pivot(int r, int s) {
        double inv = 1.0 / D[r][s];
        for (int i = 0; i < m + 2; i++) {
            if (i == r) continue;
            for (int j = 0; j < n + 2; j++) {
                if (j == s) continue;
                D[i][j] -= D[r][j] * D[i][s] * inv;
            }
        }
        for (int j = 0; j < n + 2; j++) if (j != s) D[r][j] *= inv;
        for (int i = 0; i < m + 2; i++) if (i != r) D[i][s] *= -inv;
        D[r][s] = inv;
        swap(B[r], N[s]);
    }

    bool runPhase(int phase) {
        int obj_row = (phase == 1) ? m + 1 : m;
        int iter = 0;
        cout << "\n--- Simplex Phase " << phase << " start ---\n";
        logBasis(phase, iter);

        while (true) {
            // Most negative reduced cost
            int s = -1;
            for (int j = 0; j <= n; j++) {
                if (phase == 2 && N[j] == -1) continue;
                if (s == -1 || D[obj_row][j] < D[obj_row][s] - EPS_S ||
                    (fabs(D[obj_row][j] - D[obj_row][s]) < EPS_S && N[j] < N[s]))
                    s = j;
            }

            if (D[obj_row][s] > -EPS_S) {
                cout << "  => Optimal (no negative reduced cost)\n";
                return true;
            }

            cout << "  Entering: " << varName(N[s], n)
                 << " (rc=" << fixed << setprecision(4) << D[obj_row][s] << ")\n";

            // Minimum ratio test
            int r = -1;
            double minRatio = INF_S;
            for (int i = 0; i < m; i++) {
                if (D[i][s] > EPS_S) {
                    double ratio = D[i][n + 1] / D[i][s];
                    if (ratio < minRatio - EPS_S) { minRatio = ratio; r = i; }
                }
            }
            if (r == -1) { cout << "  => Unbounded\n"; return false; }

            cout << "  Leaving:  " << varName(B[r], n)
                 << " (ratio=" << fixed << setprecision(4) << minRatio << ")\n";

            pivot(r, s);
            logBasis(phase, ++iter);

            if (iter > 10000) { cout << "  => Iteration limit (cycling?)\n"; return false; }
        }
    }

    double solve(vector<double>& solution) {
        // Check if we need phase 1
        int r = 0;
        for (int i = 1; i < m; i++)
            if (D[i][n + 1] < D[r][n + 1]) r = i;

        if (D[r][n + 1] < -EPS_S) {
            cout << "[Simplex] Negative RHS at row " << r << " — starting Phase 1.\n";
            pivot(r, n);
            if (!runPhase(1) || D[m + 1][n + 1] > EPS_S) {
                cout << "[Simplex] Infeasible.\n";
                return INF_S;
            }
        } else {
            cout << "[Simplex] Initial BFS feasible — skipping Phase 1.\n";
        }

        if (!runPhase(2)) { cout << "[Simplex] Unbounded.\n"; return INF_S; }

        solution.assign(n, 0.0);
        for (int i = 0; i < m; i++)
            if (B[i] < n) solution[B[i]] = D[i][n + 1];

        return D[m][n + 1];
    }
};

const double Simplex::EPS_S = 1e-9;
const double Simplex::INF_S = 1e18;

// ============================================================
// Big-M simplex solver for transportation LP
// Uses m+n-1 constraints (drops redundant last demand row)
// Prints verbose iteration log + verification
// ============================================================
inline double solveTransportSimplexBigM(const vector<double>& supply,
                                         const vector<double>& demand,
                                         const vector<vector<double>>& cost,
                                         vector<vector<double>>& planOut) {
    const double BIG_M_VAL = 100000.0;
    const double eps       = 1e-9;

    int m = (int)supply.size();
    int n = (int)demand.size();

    cout << "  Variables: " << m * n << " (x_ij)\n";
    cout << "  Constraints: " << m + n - 1
         << " equality (rank = m+n-1 = " << m + n - 1 << ")\n";
    cout << "  Method: Big-M with artificial variables\n\n";

    int ncon  = m + n - 1;
    int nvar  = m * n;
    int total = nvar + ncon;
    int C     = total + 1;

    // Build constraint matrix A and RHS b
    vector<vector<double>> A(ncon, vector<double>(nvar, 0));
    vector<double> b(ncon);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) A[i][i * n + j] = 1.0;
        b[i] = supply[i];
    }
    for (int j = 0; j < n - 1; j++) {
        for (int i = 0; i < m; i++) A[m + j][i * n + j] = 1.0;
        b[m + j] = demand[j];
    }

    // Cost vector: original vars + Big-M for artificials
    vector<double> cvec(total);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            cvec[i * n + j] = cost[i][j];
    for (int k = 0; k < ncon; k++)
        cvec[nvar + k] = BIG_M_VAL;

    auto varName = [&](int idx) -> string {
        if (idx < nvar) return "x" + to_string(idx / n + 1) + to_string(idx % n + 1);
        return "a" + to_string(idx - nvar + 1);
    };

    // Tableau setup: [A | I_art | b]
    vector<vector<double>> tab(ncon, vector<double>(C, 0));
    vector<int> basis(ncon);
    for (int i = 0; i < ncon; i++) {
        for (int j = 0; j < nvar; j++) tab[i][j] = A[i][j];
        tab[i][nvar + i] = 1.0;
        tab[i][C - 1]    = b[i];
        basis[i]         = nvar + i;
    }

    cout << "--- Simplex start ---\n";
    int iter = 0;

    while (true) {
        // Reduced costs: rc[j] = c[j] - c_B * col_j
        vector<double> rc(total);
        for (int j = 0; j < total; j++) {
            double s = 0;
            for (int i = 0; i < ncon; i++) s += cvec[basis[i]] * tab[i][j];
            rc[j] = cvec[j] - s;
        }

        // Real objective (artificial vars excluded)
        double real_z = 0;
        for (int i = 0; i < ncon; i++)
            if (basis[i] < nvar) real_z += cvec[basis[i]] * tab[i][C - 1];

        bool has_art = false;
        for (int i = 0; i < ncon; i++)
            if (basis[i] >= nvar && tab[i][C - 1] > eps) has_art = true;

        cout << "  [iter " << iter << "]"
             << (has_art ? " Phase 1" : " Phase 2")
             << "  obj=" << fixed << setprecision(2) << real_z << "\n";

        cout << "    BFS: ";
        int cnt = 0;
        for (int i = 0; i < ncon; i++) {
            if (tab[i][C - 1] > eps) {
                if (cnt > 0 && cnt % 6 == 0) cout << "\n         ";
                cout << varName(basis[i]) << "=" << fixed << setprecision(2)
                     << tab[i][C - 1] << "  ";
                cnt++;
            }
        }
        cout << "\n";

        // Find entering variable (most negative reduced cost)
        int entering = -1;
        double min_rc = -eps;
        for (int j = 0; j < total; j++) {
            bool is_basic = false;
            for (int k = 0; k < ncon; k++) if (basis[k] == j) { is_basic = true; break; }
            if (is_basic) continue;
            if (rc[j] < min_rc) { min_rc = rc[j]; entering = j; }
        }
        if (entering == -1) { cout << "  => Optimal (no negative reduced cost)\n"; break; }

        // Find leaving variable (minimum ratio test)
        int leaving = -1;
        double min_ratio = 1e18;
        for (int i = 0; i < ncon; i++) {
            if (tab[i][entering] > eps) {
                double ratio = tab[i][C - 1] / tab[i][entering];
                if (ratio < min_ratio - eps) { min_ratio = ratio; leaving = i; }
            }
        }
        if (leaving == -1) { cout << "  => Unbounded!\n"; return -1; }

        cout << "  Entering: " << varName(entering)
             << " (rc=" << fixed << setprecision(4) << min_rc << ")\n";
        cout << "  Leaving:  " << varName(basis[leaving])
             << " (ratio=" << fixed << setprecision(4) << min_ratio << ")\n";

        // Pivot
        double piv = tab[leaving][entering];
        for (int j = 0; j < C; j++) tab[leaving][j] /= piv;
        for (int i = 0; i < ncon; i++) {
            if (i == leaving) continue;
            double f = tab[i][entering];
            for (int j = 0; j < C; j++) tab[i][j] -= f * tab[leaving][j];
        }
        basis[leaving] = entering;
        iter++;
    }

    // Extract solution
    vector<double> x(nvar, 0);
    for (int i = 0; i < ncon; i++)
        if (basis[i] < nvar) x[basis[i]] = tab[i][C - 1];

    double total_cost = 0;
    for (int j = 0; j < nvar; j++) total_cost += cvec[j] * x[j];

    cout << "\n=== Simplex Result ===\n";
    cout << "Minimal cost: " << fixed << setprecision(2) << total_cost << "\n\n";
    cout << "Optimal plan:\n";
    planOut.assign(m, vector<double>(n));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            planOut[i][j] = x[i * n + j];
            cout << setw(8) << fixed << setprecision(2) << planOut[i][j];
        }
        cout << "\n";
    }

    cout << "\nVerification:\n  Supply: ";
    for (int i = 0; i < m; i++) {
        double s = 0;
        for (int j = 0; j < n; j++) s += x[i * n + j];
        cout << " A" << i + 1 << "=" << s << "/" << supply[i];
    }
    cout << "\n  Demand: ";
    for (int j = 0; j < n; j++) {
        double s = 0;
        for (int i = 0; i < m; i++) s += x[i * n + j];
        cout << " B" << j + 1 << "=" << s << "/" << demand[j];
    }
    cout << "\n\n";

    return total_cost;
}

// Build simplex LP from transportation data and solve
inline double solveTransportSimplex(const vector<double>& supply,
                                    const vector<double>& demand,
                                    const vector<vector<double>>& cost,
                                    vector<vector<double>>& planOut) {
    int m = (int)supply.size();
    int n = (int)demand.size();
    int vars = m * n;

    vector<vector<double>> A;
    vector<double> b;

    for (int i = 0; i < m; i++) {
        vector<double> row(vars, 0.0);
        for (int j = 0; j < n; j++) row[i * n + j] = 1.0;
        A.push_back(row);
        b.push_back(supply[i]);
    }
    for (int j = 0; j < n; j++) {
        vector<double> row(vars, 0.0);
        for (int i = 0; i < m; i++) row[i * n + j] = 1.0;
        A.push_back(row);
        b.push_back(demand[j]);
    }

    vector<double> c(vars);
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            c[i * n + j] = cost[i][j];

    Simplex solver(A, b, c);
    vector<double> x;
    double minCost = solver.solve(x);

    planOut.assign(m, vector<double>(n));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            planOut[i][j] = x[i * n + j];

    return minCost;
}
