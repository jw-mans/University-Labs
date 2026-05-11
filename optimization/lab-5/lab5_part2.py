"""
Laboratornaya rabota 5, Chast 2
Metod otsekayuschikh giperploskostey dlya nelineynoy tseli

Zadacha:
    min  f(x1, x2) = x1 + 3*x2 + exp(x1^2 + x2^2)

Ogranicheniya (dva stsenariya):

  Stsenariy A -- minimum VNUTRI dopustimoy oblasti:
    phi_1(x) = x1^2 + x2^2 - 2.0 <= 0   (krug radiusa sqrt(2) ~= 1.41)
    phi_2(x) = (x1+1)^2 + (x2+2)^2 - 9 <= 0

  Stsenariy B -- minimum NA GRANITSE dopustimoy oblasti:
    phi_1(x) = x1^2 + x2^2 - 0.25 <= 0  (krug radiusa 0.5)
    phi_2(x) = (x1+1)^2 + (x2+2)^2 - 9 <= 0

----------------------------------------------------------------------
Bezyuslovnyy minimum f:
    df/dx1 = 1 + 2*x1*E = 0  =>  x1 = -1/(2E)
    df/dx2 = 3 + 2*x2*E = 0  =>  x2 = -3/(2E)
    E = exp(x1^2+x2^2) = exp(5/(2E^2))
    Chisllenno: E ~= 1.943, x1* ~= -0.257, x2* ~= -0.772
    ||x*||^2 ~= 0.664

  Stsenariy A: ||x*||^2 = 0.664 < 2.0  =>  x* vnutri phi_1  => optim = bezuslovnyy
  Stsenariy B: ||x*||^2 = 0.664 > 0.25 =>  x* snaruzhi phi_1 => optim na granitse

----------------------------------------------------------------------
ALGORITM (Epigrafoye formulirovanie):
    min  z  s.t. f(x) <= z,  phi_i(x) <= 0
    Peremennye: (x1, x2, z), lineynaya tsel: c = (0, 0, 1)

Na kazhdoy iteratsii pri tochke x_k dobavlyaem:
  1) Otchechku dlya tseli (s parametrom alpha):
         f(x_k) + (1/alpha)*nabla_f(x_k)^T (x - x_k) <= z
     =>  [gf1/alpha, gf2/alpha, -1/alpha] * [x1,x2,z]^T
                            <= (nabla_f^T x_k - f_k) / alpha

  2) Otchechki dlya narushennykh ogranicheniy phi_i(x_k) > 0:
         phi_i(x_k) + (1/alpha)*nabla_phi_i(x_k)^T(x-x_k) <= 0
     =>  [gp1/alpha, gp2/alpha, 0] * [x1,x2,z]^T
                            <= (nabla_phi_i^T x_k - phi_i(x_k)) / alpha

Ostanov: |f(x_k) - z_k| < eps  I  vse phi_i(x_k) <= 0
----------------------------------------------------------------------
"""

import numpy as np
from scipy.optimize import linprog


def f(x):
    return x[0] + 3.0 * x[1] + np.exp(x[0]**2 + x[1]**2)


def grad_f(x):
    E = np.exp(x[0]**2 + x[1]**2)
    return np.array([1.0 + 2.0 * x[0] * E,
                     3.0 + 2.0 * x[1] * E])


def phi1(x, r2):
    return x[0]**2 + x[1]**2 - r2


def grad_phi1(x):
    return np.array([2.0 * x[0], 2.0 * x[1]])


def phi2(x):
    return (x[0] + 1.0)**2 + (x[1] + 2.0)**2 - 9.0


def grad_phi2(x):
    return np.array([2.0 * (x[0] + 1.0),
                     2.0 * (x[1] + 2.0)])


def cutting_plane_part2(r2, eps, alpha=1.0, verbose=True):
    """
    Metod otsekayuschikh giperploskostey (epigrafonaya formulirovka).

    r2    -- radius^2 dlya phi_1 (2.0 => stsenariy A, 0.25 => stsenariy B)
    eps   -- tochnost
    alpha -- koeffitsient masshtabovaniya subgradienta
    """
    c_lp = np.array([0.0, 0.0, 1.0])  # min z
    bounds = [(-3.0, 3.0), (-3.0, 3.0), (-20.0, 500.0)]

    max_iter = 300
    A_cuts = []
    b_cuts = []

    x_k = np.array([0.0, 0.0])

    def phi_all(x):
        return phi1(x, r2), phi2(x)

    def solve_lp():
        if not A_cuts:
            return linprog(c_lp, bounds=bounds, method='highs')
        return linprog(c_lp,
                       A_ub=np.array(A_cuts), b_ub=np.array(b_cuts),
                       bounds=bounds, method='highs')

    if verbose:
        print("  %4s | %9s | %9s | %10s | %10s | %8s | info"
              % ("k", "x1", "x2", "f(x)", "z", "|f-z|"))
        print("  " + "-" * 70)

    z_k = None
    x_prev = None

    for k in range(max_iter):
        f_k = f(x_k)
        gf = grad_f(x_k)

        # Otchechka dlya tseli s parametrom alpha
        a_f = np.array([gf[0] / alpha, gf[1] / alpha, -1.0 / alpha])
        b_f = (gf @ x_k - f_k) / alpha
        A_cuts.append(a_f)
        b_cuts.append(b_f)

        # Otchechki dlya narushennykh ogranicheniy
        p1, p2 = phi_all(x_k)
        gp1, gp2 = grad_phi1(x_k), grad_phi2(x_k)
        for pv, gp in [(p1, gp1), (p2, gp2)]:
            if pv > -eps * 0.1:
                a_p = np.array([gp[0] / alpha, gp[1] / alpha, 0.0])
                b_p = (gp @ x_k - pv) / alpha
                A_cuts.append(a_p)
                b_cuts.append(b_p)

        res = solve_lp()
        if not res.success:
            if verbose:
                print("  [k=%d] Oshibka LP: %s" % (k, res.message))
            break

        x_new = res.x[:2]
        z_new = res.x[2]
        f_new = f(x_new)
        gap = abs(f_new - z_new)
        p1n, p2n = phi_all(x_new)
        feasible = (p1n <= eps) and (p2n <= eps)

        if verbose:
            info = "OK" if feasible else ("phi=(%.3f,%.3f)" % (p1n, p2n))
            print("  %4d | %9.4f | %9.4f | %10.4f | %10.4f | %8.5f | %s"
                  % (k, x_new[0], x_new[1], f_new, z_new, gap, info))

        if gap < eps and feasible:
            if verbose:
                print("  [k=%d] |f - z| = %.2e < eps --> soshlos!" % (k, gap))
            x_k = x_new
            z_k = z_new
            break

        x_prev = x_k.copy()
        x_k = x_new
        z_k = z_new

        if x_prev is not None and np.linalg.norm(x_k - x_prev) < eps * 0.01:
            if verbose:
                print("  [k=%d] Maloe peremeschenie --> ostanov." % k)
            break

    return x_k, f(x_k), z_k


if __name__ == "__main__":
    sep = "=" * 65

    print(sep)
    print("  LAB. 5, CHAST 2 -- METOD OTSEKAYUSCHIKH GIPERPLOSKOSTEY")
    print(sep)
    print("  Tsel:   f(x1,x2) = x1 + 3*x2 + exp(x1^2 + x2^2)")
    print("  phi_1:  x1^2 + x2^2 - r^2 <= 0           (parametr r^2)")
    print("  phi_2:  (x1+1)^2 + (x2+2)^2 - 9 <= 0")
    print()
    print("  Bezuslovnyy minimum: x1* ~= -0.257,  x2* ~= -0.772")
    print("  ||x*||^2 ~= 0.664")
    print()

    scenarios = [
        ("A", 2.0,  "MINIMUM VNUTRI  (r^2=2.0 > 0.664)"),
        ("B", 0.25, "MINIMUM NA GRANITSE  (r^2=0.25 < 0.664)"),
    ]

    for scenario, r2, label in scenarios:
        print(sep)
        print("  STSENARIY %s: %s" % (scenario, label))
        print(sep)

        for eps in [0.1, 0.01, 0.001]:
            print("\n  --- eps = %.3f ---" % eps)
            x_opt, f_opt, z_opt = cutting_plane_part2(
                r2=r2, eps=eps, alpha=1.0, verbose=True)
            print()
            print("  REZULTAT:  x* ~= (%.4f, %.4f)" % (x_opt[0], x_opt[1]))
            print("             f* ~= %.4f" % f_opt)
            print("  Proverka ogranicheniy (dopusk = eps = %.3f):" % eps)
            v1 = phi1(x_opt, r2)
            v2 = phi2(x_opt)
            print("    phi_1(x*) = %.5f  %s" % (v1, "<= 0 OK" if v1 <= eps else ">> NARUSHENO"))
            print("    phi_2(x*) = %.5f  %s" % (v2, "<= 0 OK" if v2 <= eps else ">> NARUSHENO"))
            print()
