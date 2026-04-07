#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include <string>

#include "transport.h"
#include "simplex.h"

using namespace std;

// ============================================================
// Pretty printing helpers
// ============================================================
static void separator(const string& title) {
    cout << "\n";
    cout << "============================================================\n";
    cout << "  " << title << "\n";
    cout << "============================================================\n";
}

// ============================================================
// Shared problem data
// ============================================================
static vector<double>         SUPPLY = {8, 3, 11, 5};
static vector<double>         DEMAND = {4, 5, 4, 9, 5};
static vector<vector<double>> COST   = {
    { 8,  2,  5,  3, 14},
    {10,  4,  5,  7, 15},
    { 5,  1,  2,  1, 10},
    { 6,  3,  2,  4, 15}
};

// Task 4/5: deficit case — each supplier loses 1 unit
static vector<double> SUPPLY_DEFICIT = {7, 2, 10, 4};   // original minus 1 each
static vector<double> RATINGS        = {3, 1, 4, 2};     // supplier ratings (higher = more important)

// ============================================================
// Task 1: NWC + MODI (potentials), balanced problem
// ============================================================
static void task1_nwc_modi() {
    separator("TASK 1: Initial plan: NW-Corner | Optimisation: MODI (potentials)");

    TransportProblem prob(SUPPLY, DEMAND, COST);
    cout << "Problem data:\n";
    prob.print();
    cout << "\n";

    TransportSolver solver(prob, make_unique<NorthWestCorner>());
    solver.solve();
}

// ============================================================
// Task 2: Simplex method (Big-M), same balanced problem
// ============================================================
static void task2_simplex() {
    separator("TASK 2: Simplex Method (Big-M, Minimization)");

    cout << "Building LP from transportation data...\n";

    vector<vector<double>> plan;
    solveTransportSimplexBigM(SUPPLY, DEMAND, COST, plan);
}

// ============================================================
// Task 4: Deficit unbalanced problem (supply < demand)
//         Balanced by adding a dummy supplier with zero costs
// ============================================================
static void task4_deficit_unbalanced() {
    separator("TASK 4: Deficit Unbalanced Problem (dummy supplier)");

    TransportProblem deficitProb(SUPPLY_DEFICIT, DEMAND, COST);

    cout << "Unbalanced problem data:\n";
    deficitProb.print();
    cout << "\n";

    cout << "Supply=" << deficitProb.totalSupply()
         << " < Demand=" << deficitProb.totalDemand()
         << " — adding dummy supplier to cover deficit of "
         << deficitProb.totalDemand() - deficitProb.totalSupply() << "\n\n";

    TransportProblem balanced = makeDeficitBalanced(deficitProb);

    cout << "Balanced problem data:\n";
    balanced.print();
    cout << "\n";

    TransportSolver solver(balanced, make_unique<NorthWestCorner>());
    solver.solve();
}

// ============================================================
// Task 5: Priority-based assignment (Supplier Ratings)
//         Highest-rated supplier fills cheapest consumers first
// ============================================================
static void task5_rated_priority() {
    separator("TASK 5: Priority-based Distribution (Supplier Ratings)");

    TransportProblem prob(SUPPLY, DEMAND, COST);
    cout << "Problem data:\n";
    prob.print();
    cout << "\n(Higher rating = higher priority, fills cheapest consumers first)\n\n";

    GreedyPriorityTransportSolver solver(prob, RATINGS);
    solver.solve();
}

// ============================================================
// Main pipeline
// ============================================================
int main() {
    task1_nwc_modi();
    task2_simplex();
    task4_deficit_unbalanced();
    task5_rated_priority();

    return 0;
}
