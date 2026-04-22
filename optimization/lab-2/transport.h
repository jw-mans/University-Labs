#pragma once
#include <iostream>
#include <vector>
#include <numeric>
#include <limits>
#include <algorithm>
#include <set>
#include <cmath>
#include <utility>
#include <memory>
#include <stdexcept>
#include <string>

using namespace std;

const double EPS = 1e-10;

// ============================================================
// Data class
// ============================================================
struct TransportProblem {
    vector<double> supply;
    vector<double> demand;
    vector<vector<double>> cost;

    TransportProblem(vector<double> a, vector<double> b, vector<vector<double>> c)
        : supply(move(a)), demand(move(b)), cost(move(c)) {}

    int numSuppliers() const { return (int)supply.size(); }
    int numConsumers()  const { return (int)demand.size(); }

    double totalSupply() const { double s = 0; for (double v : supply) s += v; return s; }
    double totalDemand() const { double s = 0; for (double v : demand)  s += v; return s; }
    bool   isBalanced()  const { return fabs(totalSupply() - totalDemand()) < EPS; }

    void print() const {
        int m = numSuppliers(), n = numConsumers();
        cout << "  Costs:\n";
        for (int i = 0; i < m; ++i) {
            cout << "    ";
            for (int j = 0; j < n; ++j)
                cout << cost[i][j] << "\t";
            cout << "| supply=" << supply[i] << "\n";
        }
        cout << "    demand: ";
        for (double d : demand) cout << d << "\t";
        cout << "\n";
        cout << "  Total supply=" << totalSupply()
             << "  Total demand=" << totalDemand() << "\n";
    }
};

// Factory: add a dummy supplier to cover deficit (supply < demand)
inline TransportProblem makeDeficitBalanced(const TransportProblem& p) {
    double deficit = p.totalDemand() - p.totalSupply();
    if (deficit <= EPS) return p;

    vector<double> newSupply = p.supply;
    newSupply.push_back(deficit);

    vector<vector<double>> newCost = p.cost;
    newCost.push_back(vector<double>(p.numConsumers(), 0.0));

    return TransportProblem(newSupply, p.demand, newCost);
}

// ============================================================
// Strategy: Initial Plan Builder
// ============================================================
class InitialPlanBuilder {
public:
    virtual ~InitialPlanBuilder() = default;
    virtual string name() const = 0;
    virtual void build(const TransportProblem& prob,
                       vector<vector<double>>& x,
                       vector<pair<int,int>>& basis) const = 0;
};

class NorthWestCorner : public InitialPlanBuilder {
public:
    string name() const override { return "north-west corner method"; }

    void build(const TransportProblem& prob,
               vector<vector<double>>& x,
               vector<pair<int,int>>& basis) const override {
        int m = prob.numSuppliers(), n = prob.numConsumers();

        // Инициализируем план нулями
        x.assign(m, vector<double>(n, 0.0));
        basis.clear();

        // Копируем запасы и потребности — будем их уменьшать по мере заполнения
        vector<double> supply = prob.supply;
        vector<double> demand = prob.demand;

        // Начинаем с северо-западного угла (i=0, j=0)
        int i = 0, j = 0;
        while (i < m && j < n) {
            // Пропускаем поставщика, если его запас исчерпан
            if (supply[i] < EPS) { ++i; continue; }
            // Пропускаем потребителя, если его потребность удовлетворена
            if (demand[j] < EPS) { ++j; continue; }

            // Назначаем максимально возможное количество — минимум из остатка запаса и потребности
            double amount = min(supply[i], demand[j]);
            x[i][j] = amount;
            basis.push_back({i, j});  // Ячейка входит в базис
            supply[i] -= amount;
            demand[j] -= amount;

            if (supply[i] < EPS && demand[j] < EPS) {
                // Вырожденность: оба исчерпаны одновременно — добавляем нулевую базисную ячейку,
                // чтобы сохранить размер базиса m+n-1, и двигаемся по строке или столбцу
                if (i + 1 < m)      { basis.push_back({i + 1, j}); ++i; }
                else if (j + 1 < n) { basis.push_back({i, j + 1}); ++j; }
            } else if (supply[i] < EPS) {
                ++i;  // Запас поставщика i исчерпан — переходим к следующему поставщику
            } else {
                ++j;  // Потребность потребителя j удовлетворена — переходим к следующему потребителю
            }
        }

        if ((int)basis.size() != m + n - 1)
            cerr << "Warning: basis size = " << basis.size()
                 << ", required " << m + n - 1 << ". Possible issues.\n";
    }
};

// ============================================================
// Solver: Potentials Method (MODI)
// ============================================================
class TransportSolver {
protected:
    TransportProblem prob;
    vector<vector<double>> x;
    vector<pair<int,int>> basis;
    unique_ptr<InitialPlanBuilder> strategy;

    // MODI, шаг 1: вычисляем потенциалы u[i] и v[j] из условия u[i]+v[j]=c[i][j] для базисных ячеек.
    // Фиксируем u[0]=0, остальные потенциалы находим итеративно по уже известным соседям.
    bool computePotentials(vector<double>& u, vector<double>& v) const {
        int m = prob.numSuppliers(), n = prob.numConsumers();
        u.assign(m, numeric_limits<double>::quiet_NaN());
        v.assign(n, numeric_limits<double>::quiet_NaN());
        u[0] = 0.0;  // Опорное значение — один потенциал задаём свободно

        bool changed;
        do {
            changed = false;
            for (const auto& [r, c] : basis) {
                // Если известен u[r], находим v[c] = c[r][c] - u[r]
                if (!isnan(u[r]) && isnan(v[c]))  { v[c] = prob.cost[r][c] - u[r]; changed = true; }
                // Если известен v[c], находим u[r] = c[r][c] - v[c]
                else if (!isnan(v[c]) && isnan(u[r])) { u[r] = prob.cost[r][c] - v[c]; changed = true; }
            }
        } while (changed);  // Повторяем, пока хоть один потенциал не определён

        // Если остались NaN — базис несвязен, метод неприменим
        for (double ui : u) if (isnan(ui)) return false;
        for (double vj : v) if (isnan(vj)) return false;
        return true;
    }

    // MODI, шаг 2: ищем небазисную ячейку с наибольшим отрицательным оценочным числом
    // delta[i][j] = c[i][j] - (u[i] + v[j]).
    // Если все delta >= 0 — план оптимален. Иначе — берём ячейку с минимальным delta (наибольшее улучшение).
    pair<int,int> findEnteringCell(const vector<double>& u, const vector<double>& v,
                                   double& minDelta) const {
        set<pair<int,int>> basisSet(basis.begin(), basis.end());
        minDelta = 0.0;
        pair<int,int> entering = {-1, -1};

        for (int i = 0; i < prob.numSuppliers(); ++i)
            for (int j = 0; j < prob.numConsumers(); ++j) {
                if (basisSet.count({i, j})) continue;  // Пропускаем базисные ячейки
                double delta = prob.cost[i][j] - (u[i] + v[j]);  // Оценочное число
                if (delta < -EPS && (entering.first == -1 || delta < minDelta)) {
                    minDelta = delta;
                    entering = {i, j};  // Запоминаем «лучшую» входящую ячейку
                }
            }
        return entering;
    }

    // MODI, шаг 3: DFS для построения цикла пересчёта.
    // Цикл — замкнутый путь по базисным ячейкам с чередованием горизонтальных и вертикальных ходов,
    // начинающийся и заканчивающийся во входящей ячейке.
    // lastDir: -1=начало, 0=шли горизонтально, 1=шли вертикально (направление чередуется).
    static bool cycleDFS(pair<int,int> cur, pair<int,int> start, int lastDir,
                         vector<pair<int,int>>& path, const set<pair<int,int>>& avail) {
        if (cur == start && !path.empty()) return true;  // Цикл замкнулся

        // Горизонтальный ход (dir=0): ищем ячейку в той же строке
        if (lastDir != 0)
            for (const auto& cell : avail) {
                if (cell.first != cur.first || cell.second == cur.second) continue;
                if (cell == start && !path.empty()) { path.push_back(cell); return true; }
                if (find(path.begin(), path.end(), cell) == path.end()) {
                    path.push_back(cell);
                    if (cycleDFS(cell, start, 0, path, avail)) return true;
                    path.pop_back();
                }
            }

        // Вертикальный ход (dir=1): ищем ячейку в том же столбце
        if (lastDir != 1)
            for (const auto& cell : avail) {
                if (cell.second != cur.second || cell.first == cur.first) continue;
                if (cell == start && !path.empty()) { path.push_back(cell); return true; }
                if (find(path.begin(), path.end(), cell) == path.end()) {
                    path.push_back(cell);
                    if (cycleDFS(cell, start, 1, path, avail)) return true;
                    path.pop_back();
                }
            }

        return false;
    }

    // MODI, шаг 4: запускаем DFS, формируем полный цикл начиная с входящей ячейки.
    vector<pair<int,int>> buildCycle(pair<int,int> entering) const {
        // Доступные узлы цикла — все базисные ячейки плюс входящая
        set<pair<int,int>> avail(basis.begin(), basis.end());
        avail.insert(entering);

        vector<pair<int,int>> path;
        if (!cycleDFS(entering, entering, -1, path, avail)) {
            cerr << "Error: failed to build recalculation cycle!\n";
            return {};
        }

        // Собираем цикл: входящая ячейка + найденный путь (без дублирования конечной точки)
        vector<pair<int,int>> cycle = {entering};
        cycle.insert(cycle.end(), path.begin(), path.end());
        if (cycle.size() > 1 && cycle.back() == entering) cycle.pop_back();
        return cycle;
    }

    // MODI, шаг 5: находим theta — максимальное количество, которое можно перераспределить.
    // theta = минимум из значений x в "минусовых" ячейках цикла (нечётные позиции).
    double findTheta(const vector<pair<int,int>>& cycle) const {
        double theta = numeric_limits<double>::max();
        for (size_t k = 1; k < cycle.size(); k += 2)
            theta = min(theta, x[cycle[k].first][cycle[k].second]);
        return theta;
    }

    // MODI, шаг 6: перераспределяем поставки вдоль цикла на theta.
    // Чётные ячейки цикла ("+") получают +theta, нечётные ("-") теряют -theta.
    // Ячейка, обнулившаяся в минусах, выходит из базиса; входящая — входит.
    void redistribute(const vector<pair<int,int>>& cycle, double theta, pair<int,int> entering) {
        for (size_t k = 0; k < cycle.size(); ++k) {
            auto [i, j] = cycle[k];
            x[i][j] += (k % 2 == 0) ? theta : -theta;
        }

        // Ищем выходящую ячейку — ту, которая стала нулевой в минусах (нечётные позиции цикла)
        pair<int,int> leaving = {-1, -1};
        for (size_t k = 1; k < cycle.size(); k += 2)
            if (x[cycle[k].first][cycle[k].second] < EPS) { leaving = cycle[k]; break; }

        if (leaving.first == -1) { cerr << "Error: leaving cell not found (theta=" << theta << ")\n"; return; }

        // Удаляем выходящую ячейку из базиса и добавляем входящую
        auto it = find(basis.begin(), basis.end(), leaving);
        if (it != basis.end()) basis.erase(it);
        else cerr << "Warning: leaving cell not found in basis\n";

        basis.push_back(entering);
    }

    double totalCost() const {
        double cost = 0.0;
        for (int i = 0; i < prob.numSuppliers(); ++i)
            for (int j = 0; j < prob.numConsumers(); ++j)
                cost += x[i][j] * prob.cost[i][j];
        return cost;
    }

    void printPlan() const {
        for (const auto& row : x) {
            for (double val : row) cout << val << "\t";
            cout << "\n";
        }
    }

public:
    TransportSolver(TransportProblem p, unique_ptr<InitialPlanBuilder> s)
        : prob(move(p)), strategy(move(s)) {}

    virtual ~TransportSolver() = default;

    virtual void solve() {
        // Метод MODI (метод потенциалов) работает только с балансированной задачей
        if (!prob.isBalanced()) {
            cerr << "Error: problem is not balanced (supply=" << prob.totalSupply()
                 << ", demand=" << prob.totalDemand() << ")\n";
            return;
        }

        // Шаг 0: строим начальный допустимый план выбранным методом (например, СЗУ)
        strategy->build(prob, x, basis);
        cout << "Initial plan (" << strategy->name() << "):\n";
        printPlan();
        cout << "Initial cost: " << totalCost() << "\n\n";

        int iter = 0;
        while (true) {
            // Шаг 1: вычисляем потенциалы u[i], v[j] по текущему базису
            vector<double> u, v;
            if (!computePotentials(u, v)) {
                cerr << "Error: failed to compute potentials (basis disconnected).\n";
                break;
            }

            // Шаг 2: проверяем критерий оптимальности — ищем ячейку с delta < 0
            double minDelta;
            auto entering = findEnteringCell(u, v, minDelta);
            if (entering.first == -1) { 
                cout << "Optimal plan reached!\n"; 
                break; 
            }  // Все delta >= 0 — оптимум

            cout << "Iter " << ++iter << ": entering (" << entering.first << ","
                 << entering.second << ") delta=" << minDelta << "\n";

            // Шаг 3-4: строим цикл пересчёта через входящую ячейку
            auto cycle = buildCycle(entering);
            if (cycle.empty()) break;

            cout << "  Cycle: ";
            for (size_t k = 0; k < cycle.size(); ++k) {
                cout << "(" << cycle[k].first << "," << cycle[k].second << ")";
                if (k + 1 < cycle.size()) cout << " -> ";
            }
            cout << " -> (" << cycle[0].first << "," << cycle[0].second << ")\n";

            // Шаг 5: находим величину перераспределения
            double theta = findTheta(cycle);
            cout << "  theta=" << theta << "\n";

            // Шаг 6: перераспределяем поставки, обновляем базис
            redistribute(cycle, theta, entering);

            cout << "  Plan after redistribution:\n";
            printPlan();
            cout << "  Cost: " << totalCost() << "\n\n";
        }

        cout << "-------------------------\n";
        cout << "Optimal plan:\n";
        printPlan();
        cout << "Total cost: " << totalCost() << "\n";
    }
};

// ============================================================
// Derived Solver: surplus unbalanced problem with supplier ratings
// Higher rating = higher priority (cargo taken from lowest-rated first)
// ============================================================
class RatedTransportSolver : public TransportSolver {
    static TransportProblem balance(const TransportProblem& p, const vector<double>& ratings) {
        if ((int)ratings.size() != p.numSuppliers())
            throw invalid_argument("Error: ratings size != number of suppliers.");
        if (p.totalSupply() < p.totalDemand() - EPS)
            throw runtime_error("Error: supply < demand. No solution.");
        if (fabs(p.totalSupply() - p.totalDemand()) < EPS) {
            cout << "Already balanced.\n";
            return p;
        }

        double surplus = p.totalSupply() - p.totalDemand();
        cout << "Surplus to remove: " << surplus << "\n";

        // Sort suppliers ascending by rating (lowest rated = cut first)
        vector<int> idx(p.numSuppliers());
        for (int i = 0; i < p.numSuppliers(); ++i) idx[i] = i;
        sort(idx.begin(), idx.end(), [&](int a, int b) { return ratings[a] < ratings[b]; });

        vector<double> balanced = p.supply;
        for (int i : idx) {
            if (surplus <= 0) break;
            double take = min(balanced[i], surplus);
            balanced[i] -= take;
            surplus -= take;
            cout << "  Supplier " << i + 1 << " (rating=" << ratings[i]
                 << "): cut " << take << ", remaining=" << balanced[i] << "\n";
        }
        if (surplus > EPS)
            throw runtime_error("Error: could not eliminate surplus.");

        // Remove suppliers with zero supply — they cause disconnected basis
        vector<double> filteredSupply;
        vector<vector<double>> filteredCost;
        for (int i = 0; i < p.numSuppliers(); ++i) {
            if (balanced[i] > EPS) {
                filteredSupply.push_back(balanced[i]);
                filteredCost.push_back(p.cost[i]);
            } else {
                cout << "  Supplier " << i + 1 << " fully cut — removed from problem.\n";
            }
        }

        cout << "Balanced supply: ";
        for (double v : filteredSupply) cout << v << " ";
        cout << "\n\n";

        return TransportProblem(filteredSupply, p.demand, filteredCost);
    }

public:
    RatedTransportSolver(const TransportProblem& p, const vector<double>& ratings,
                         unique_ptr<InitialPlanBuilder> s)
        : TransportSolver(balance(p, ratings), move(s)) {}
};

// ============================================================
// Greedy Priority Solver: assigns cargo from highest-rated
// supplier first; each supplier fills cheapest consumers first.
// Does not optimise — produces a priority-feasible plan.
// ============================================================
class GreedyPriorityTransportSolver {
    TransportProblem prob;
    vector<double>   ratings;

public:
    GreedyPriorityTransportSolver(TransportProblem p, vector<double> r)
        : prob(move(p)), ratings(move(r)) {}

    void solve() {
        int m = prob.numSuppliers(), n = prob.numConsumers();

        // Sort suppliers descending by rating
        vector<int> order(m);
        iota(order.begin(), order.end(), 0);
        sort(order.begin(), order.end(), [&](int a, int b) {
            return ratings[a] != ratings[b] ? ratings[a] > ratings[b] : a < b;
        });

        cout << "Ratings: ";
        for (int i = 0; i < m; i++) cout << "A" << i + 1 << "=" << ratings[i] << "  ";
        cout << "\nPriority order (highest first): ";
        for (int k = 0; k < m; k++) {
            cout << "A" << order[k] + 1 << "(r=" << ratings[order[k]] << ")";
            if (k + 1 < m) cout << " > ";
        }
        cout << "\n\n";

        vector<vector<double>> plan(m, vector<double>(n, 0));
        vector<double> rem_demand = prob.demand;
        vector<double> rem_supply = prob.supply;

        cout << "--- Priority distribution ---\n\n";

        int stage = 0;
        for (int i : order) {
            ++stage;
            cout << "Stage " << stage << ": A" << i + 1
                 << " (supply=" << rem_supply[i]
                 << ", rating=" << ratings[i] << ")\n";

            // Sort consumers by cost for this supplier
            vector<int> cons_order(n);
            iota(cons_order.begin(), cons_order.end(), 0);
            sort(cons_order.begin(), cons_order.end(), [&](int j1, int j2) {
                return prob.cost[i][j1] != prob.cost[i][j2]
                    ? prob.cost[i][j1] < prob.cost[i][j2]
                    : j1 < j2;
            });

            cout << "  Cost order for A" << i + 1 << ": ";
            for (int j : cons_order)
                cout << "B" << j + 1 << "(c=" << prob.cost[i][j] << ") ";
            cout << "\n";

            bool any = false;
            for (int j : cons_order) {
                if (rem_supply[i] < EPS) break;
                if (rem_demand[j] < EPS) {
                    cout << "  B" << j + 1 << ": demand exhausted, skip\n";
                    continue;
                }
                double amount = min(rem_supply[i], rem_demand[j]);
                plan[i][j]    = amount;
                rem_supply[i] -= amount;
                rem_demand[j] -= amount;
                any = true;
                cout << "  -> x(" << i + 1 << "," << j + 1 << ") = " << amount
                     << "  |  A" << i + 1 << " remaining: " << rem_supply[i]
                     << ",  B" << j + 1 << " remaining: " << rem_demand[j] << "\n";
            }
            if (!any) cout << "  (no assignments made)\n";

            cout << "  Remaining demand: (";
            for (int j = 0; j < n; j++) {
                cout << rem_demand[j];
                if (j + 1 < n) cout << ", ";
            }
            cout << ")\n\n";
        }

        // Print result
        cout << "Priority assignment plan:\n";
        double total_cost = 0;
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout << setw(6) << fixed << setprecision(0) << plan[i][j];
                total_cost += plan[i][j] * prob.cost[i][j];
            }
            cout << "\n";
        }
        cout << "Total cost: " << fixed << setprecision(2) << total_cost << "\n\n";

        cout << "Cost breakdown:\n";
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                if (plan[i][j] > EPS)
                    cout << "  x(" << i + 1 << "," << j + 1 << ")="
                         << fixed << setprecision(0) << plan[i][j]
                         << " * c=" << prob.cost[i][j]
                         << " = " << plan[i][j] * prob.cost[i][j] << "\n";
        cout << "  Sum = " << fixed << setprecision(2) << total_cost << "\n";
    }
};
