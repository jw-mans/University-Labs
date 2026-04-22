/*
 * Lab 3 (OOP refactor)
 * 1D minimization: bisection (dichotomy) and golden section
 *
 * f(x) = e^x - x^3/3 + 2x   on  [-2, -1]
 *
 * Design patterns applied:
 *   Strategy       -- Minimizer1D / BisectionMinimizer / GoldenSectionMinimizer
 *   Template Method -- Minimizer1D::run() orchestrates, solve() is overridden
 *   Value Object   -- MinResult
 *   Factory        -- make_minimizer()
 */

#include <cmath>
#include <cstdio>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// ─── Function wrapper with evaluation counter ─────────────────────────────────

class Function1D {
public:
    using Fn = std::function<double(double)>;

    Function1D(Fn f, Fn df = {}, Fn d2f = {})
        : f_(std::move(f)), df_(std::move(df)), d2f_(std::move(d2f)) {}

    double operator()(double x)   { ++calls_; return f_(x); }  // counted
    double eval(double x)   const { return f_(x); }            // uncounted
    double deriv(double x)  const { return df_(x); }
    double deriv2(double x) const { return d2f_(x); }

    int  calls() const noexcept { return calls_; }
    void reset()       noexcept { calls_ = 0; }

private:
    Fn  f_, df_, d2f_;
    int calls_ = 0;
};

// ─── Result value object ──────────────────────────────────────────────────────

struct MinResult {
    double x_min;
    double f_min;
    int    iterations;
    int    calls;
    int    analytic;
};

// ─── Abstract base: Strategy interface + Template Method ─────────────────────

class Minimizer1D {
public:
    virtual ~Minimizer1D() = default;

    // Template Method: resets counter, runs algorithm, assembles result
    MinResult run(Function1D& fn, double a, double b, double eps) const {
        fn.reset();
        auto [x, iter] = solve(fn, a, b, eps);
        return { x, fn.eval(x), iter, fn.calls(), analytic_calls(a, b, eps) };
    }

    virtual std::string name()                                  const = 0;
    virtual int analytic_calls(double a, double b, double eps)  const = 0;

protected:
    // Returns {x_min, iteration_count}
    virtual std::pair<double, int> solve(
        Function1D& fn, double a, double b, double eps)         const = 0;
};

// ─── Bisection (dichotomy) strategy ──────────────────────────────────────────

class BisectionMinimizer final : public Minimizer1D {
public:
    explicit BisectionMinimizer(double delta_ratio = 0.1)
        : delta_ratio_(delta_ratio) {}

    std::string name() const override { return "Bisection (delta=eps/10)"; }

    int analytic_calls(double a, double b, double eps) const override {
        return 2 * static_cast<int>(std::ceil(std::log2((b - a) / eps)));
    }

protected:
    std::pair<double, int> solve(
        Function1D& fn, double a, double b, double eps) const override
    {
        int iter = 0;
        const double delta = eps * delta_ratio_;
        while (b - a > eps) {
            double mid = (a + b) / 2.0;
            if (fn(mid - delta) < fn(mid + delta)) b = mid + delta;
            else                                    a = mid - delta;
            ++iter;
        }
        return { (a + b) / 2.0, iter };
    }

private:
    double delta_ratio_;
};

// ─── Golden section strategy ──────────────────────────────────────────────────

class GoldenSectionMinimizer final : public Minimizer1D {
public:
    std::string name() const override { return "Golden section"; }

    int analytic_calls(double a, double b, double eps) const override {
        static const double PHI = (1.0 + std::sqrt(5.0)) / 2.0;
        return static_cast<int>(std::ceil(std::log((b - a) / eps) / std::log(PHI))) + 1;
    }

protected:
    std::pair<double, int> solve(
        Function1D& fn, double a, double b, double eps) const override
    {
        static const double PHI = (1.0 + std::sqrt(5.0)) / 2.0;
        static const double TAU = PHI - 1.0;        // ≈ 0.618
        static const double INV = 1.0 - TAU;        // ≈ 0.382

        double x1 = a + INV * (b - a), f1 = fn(x1);
        double x2 = a + TAU * (b - a), f2 = fn(x2);
        int iter = 0;

        while (b - a > eps) {
            if (f1 < f2) {
                b  = x2; x2 = x1; f2 = f1;
                x1 = a + INV * (b - a); f1 = fn(x1);
            } else {
                a  = x1; x1 = x2; f1 = f2;
                x2 = a + TAU * (b - a); f2 = fn(x2);
            }
            ++iter;
        }
        return { (a + b) / 2.0, iter };
    }
};

// ─── Factory ──────────────────────────────────────────────────────────────────

std::unique_ptr<Minimizer1D> make_minimizer(const std::string& type) {
    if (type == "bisection") return std::make_unique<BisectionMinimizer>();
    if (type == "golden")    return std::make_unique<GoldenSectionMinimizer>();
    throw std::invalid_argument("Unknown minimizer type: " + type);
}

// ─── Unimodality checker ──────────────────────────────────────────────────────

class UnimodalityChecker {
public:
    explicit UnimodalityChecker(const Function1D& fn, int steps = 10)
        : fn_(fn), steps_(steps) {}

    void print(double a, double b) const {
        printf("=== Unimodality: f(x)=e^x - x^3/3 + 2x on [%g, %g] ===\n", a, b);
        printf("  f'(x)  = e^x - x^2 + 2\n");
        printf("  f''(x) = e^x - 2x\n\n");
        printf("  f'(%g) = %+.4f  (< 0 => decreasing)\n", a, fn_.deriv(a));
        printf("  f'(%g) = %+.4f  (> 0 => increasing)\n", b, fn_.deriv(b));
        printf("  => f' changes sign => unique minimum exists\n\n");

        bool all_pos = true;
        double h = (b - a) / steps_;
        printf("  f''(x) on [%g, %g]:\n", a, b);
        for (int i = 0; i <= steps_; ++i) {
            double x = a + i * h;
            double v = fn_.deriv2(x);
            if (v <= 0) all_pos = false;
            printf("    f''(%+.4f) = %.4f\n", x, v);
        }
        printf("  => f''(x) %s => strictly convex => unimodal\n\n",
               all_pos ? "> 0 everywhere" : "NOT > 0 everywhere (!)");
    }

private:
    const Function1D& fn_;
    int               steps_;
};

// ─── Table printer ────────────────────────────────────────────────────────────

class TablePrinter {
public:
    void print_method(const Minimizer1D& m, Function1D& fn,
                      const std::vector<double>& epsilons, double a, double b) const
    {
        printf("=== %s ===\n", m.name().c_str());
        printf("%-8s  %-11s  %-11s  %-6s  %-6s  %-6s\n",
               "eps", "x*", "f(x*)", "iter", "calls", "N_anal");
        printf("-----------------------------------------------------------\n");
        for (double eps : epsilons) {
            MinResult r = m.run(fn, a, b, eps);
            printf("%.3f     %+.6f   %+.6f   %-6d  %-6d  %-6d\n",
                   eps, r.x_min, r.f_min, r.iterations, r.calls, r.analytic);
        }
        printf("\n");
    }

    void print_comparison(const Minimizer1D& m1, const Minimizer1D& m2,
                          Function1D& fn, const std::vector<double>& epsilons,
                          double a, double b) const
    {
        printf("=== Comparison: function evaluations ===\n");
        printf("%-8s  %-14s  %-14s  %-8s\n",
               "eps", "Bisection", "Golden sect.", "Speedup");
        printf("--------------------------------------------\n");
        for (double eps : epsilons) {
            int c1 = m1.run(fn, a, b, eps).calls;
            int c2 = m2.run(fn, a, b, eps).calls;
            printf("%.3f     %-14d  %-14d  x%.2f\n", eps, c1, c2, (double)c1 / c2);
        }
        static const double PHI = (1.0 + std::sqrt(5.0)) / 2.0;
        printf("\nAsymptotics: N_dich/N_gold -> 2*ln(phi)/ln(2) = %.4f\n",
               2.0 * std::log(PHI) / std::log(2.0));
    }

    void print_analytic(const Minimizer1D& m1, const Minimizer1D& m2,
                        const std::vector<double>& epsilons, double a, double b) const
    {
        printf("\n=== Analytical call-count formulas ===\n");
        printf("Bisection:    N = 2 * ceil( log2((b-a)/eps) )\n");
        printf("Golden sect.: N = ceil( ln((b-a)/eps) / ln(phi) ) + 1\n");
        printf("\n%-8s  %-12s  %-12s\n", "eps", "N_dich", "N_gold");
        for (double eps : epsilons)
            printf("%.3f     %-12d  %-12d\n", eps,
                   m1.analytic_calls(a, b, eps), m2.analytic_calls(a, b, eps));
    }
};

// ─── main ─────────────────────────────────────────────────────────────────────

int main() {
    const double A = -2.0, B = -1.0;
    const std::vector<double> epsilons = { 0.1, 0.01, 0.001 };

    Function1D fn(
        [](double x) { return std::exp(x) - x*x*x/3.0 + 2.0*x; },
        [](double x) { return std::exp(x) - x*x + 2.0; },
        [](double x) { return std::exp(x) - 2.0*x; }
    );

    UnimodalityChecker(fn).print(A, B);

    auto bisect = make_minimizer("bisection");
    auto golden = make_minimizer("golden");
    TablePrinter printer;

    printer.print_method(*bisect, fn, epsilons, A, B);
    printer.print_method(*golden, fn, epsilons, A, B);
    printer.print_comparison(*bisect, *golden, fn, epsilons, A, B);
    printer.print_analytic(*bisect, *golden, epsilons, A, B);

    return 0;
}
