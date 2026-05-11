"""
Laboratornaya rabota 5, Chast 1
Metod otcekayuschikh giperploskostey

Zadacha:
    min  f(x) = x1 + 2*x2 + x3 + 2*x4
    s.t. phi_1(x) = x1^2  + x2^2  + x3^2  + x4^2  - 10 <= 0  (kvadratichnoe)
         phi_2(x) = x1^4  + x2^4  + x3^4  + x4^4  - 40 <= 0  (stepen 4)
         phi_3(x) = 1.2*x1^4 + x2^4 + 1.1*x3^4 + x4^4 - 42 <= 0 (stepen 4)

Vse peremennye prisutstvuyut v kazhdom ogranichenii.
Tolko phi_1 -- kvadratichnaya funktsiya.

Analiticheskoe reshenie:
    Na phi_1 (sfera radiusa sqrt(10)) min dostigaetsya v napravlenii -c:
    x* = (-1, -2, -1, -2),  f* = -10
    Proverka: phi_1(x*)=0 (aktivno), phi_2(x*)=34-40<0, phi_3(x*)=34.3-42<0

----------------------------------------------------------------------
DOKAZATELSTVO FORMULY SUBGRADIENTA S PROIZVOLNYM alpha
----------------------------------------------------------------------
Pust phi_i: R^n -> R -- vypuklaya funktsiya, x_k -- tekuschaya tochka,
g = grad phi_i(x_k) -- subgradient (= gradient dlya gladkikh funktsiy).

Opredelenie subgradienta (vypuklaya funktsiya):
    phi_i(x) >= phi_i(x_k) + g^T (x - x_k)   dlya vsekh x         (1)

Vvedem proizvolnyy parametr alpha > 0. Iz (1) dlya dopustimoy tochki x:
    phi_i(x) <= 0  =>  phi_i(x_k) + g^T(x - x_k) <= 0
    =>  g^T x  <=  g^T x_k - phi_i(x_k)                          (2)

Delim (2) na alpha > 0 (znak neravenstva sokhranjaetsya):
    (g/alpha)^T x  <=  (g^T x_k - phi_i(x_k)) / alpha            (3)

Eto i est otsekayuschaya giperploskost s masshtabированным subgradientom.

Svoystva:
 * Ni odna dopustimaya tochka ne isklyuchaetsya: (3) vypolneno dlya vsekh
   x s phi_i(x) <= 0 (poluchen iz neravenstva subgradienta).
 * x_k otsechaetsya: podstavim x = x_k v (3):
     (g/alpha)^T x_k = g^T x_k / alpha  >  (g^T x_k - phi_i(x_k)) / alpha
   (tak kak phi_i(x_k) > 0 => pri delenii na alpha raznost menshe)
 * alpha = 1: standartnyy shag metoda otsekayuschikh giperploskostey Kelli.

Vyvod: otsekayuschaya giperploskost s parametrom alpha > 0:
    a = g / alpha,    b = (g^T x_k - phi_i(x_k)) / alpha
    a^T x <= b
korektno suzhaet dopustimoe mnozhestvo.
----------------------------------------------------------------------
"""

import numpy as np
from scipy.optimize import linprog

# ---------- Zadacha ----------
c_obj = np.array([1.0, 2.0, 1.0, 2.0])
N = 4
M = 5.0  # nachalnoe kompaktnoe mnozhestvo [-M, M]^N


def phi(x, i):
    if i == 1:
        return x[0]**2 + x[1]**2 + x[2]**2 + x[3]**2 - 10.0
    elif i == 2:
        return x[0]**4 + x[1]**4 + x[2]**4 + x[3]**4 - 40.0
    elif i == 3:
        return 1.2*x[0]**4 + x[1]**4 + 1.1*x[2]**4 + x[3]**4 - 42.0


def grad_phi(x, i):
    if i == 1:
        return np.array([2*x[0], 2*x[1], 2*x[2], 2*x[3]])
    elif i == 2:
        return np.array([4*x[0]**3, 4*x[1]**3, 4*x[2]**3, 4*x[3]**3])
    elif i == 3:
        return np.array([4.8*x[0]**3, 4*x[1]**3, 4.4*x[2]**3, 4*x[3]**3])


def cutting_plane(eps, alpha=1.0, verbose=True):
    """
    Metod otsekayuschikh giperploskostey (Kelley, 1960).

    eps   -- tochnost ostanova
    alpha -- koeffitsient masshtabovaniya subgradienta (alpha > 0)
    verbose -- pechatat detali iteratsiy

    Vozvraschaet (x_opt, f_opt, dual_vars, A_cuts, b_cuts)
    """
    m = 3
    max_iter = 500

    # Nachalnye otchechi: kub [-M, M]^N
    A_cuts = []
    b_cuts = []
    for j in range(N):
        e = np.zeros(N)
        e[j] = 1.0
        A_cuts.append(e.copy());  b_cuts.append(M)
        A_cuts.append(-e.copy()); b_cuts.append(M)

    def solve_lp():
        A = np.array(A_cuts)
        b = np.array(b_cuts)
        return linprog(c_obj, A_ub=A, b_ub=b,
                       bounds=[(None, None)] * N, method='highs')

    res = solve_lp()
    if not res.success:
        raise RuntimeError("Ne udalos reshit nachalnuyu LP.")

    x_k = res.x.copy()
    x_prev = None

    if verbose:
        print("  x_0 = %s,  f = %.5f" % (np.round(x_k, 4), c_obj @ x_k))

    for k in range(max_iter):
        phis = [phi(x_k, i + 1) for i in range(m)]
        max_phi = max(phis)

        # Kriteriy 1: x_k dopustimo
        if max_phi <= 1e-8:
            if verbose:
                print("  [k=%d] x_k in Omega --> reshenie naydeno!" % k)
            break

        # Kriteriy 2: maloe peremeschenie
        if x_prev is not None:
            diff = np.linalg.norm(x_k - x_prev)
            if verbose:
                print("  [k=%d] ||dx|| = %.6f,  phi_max = %.6f" % (k, diff, max_phi))
            if diff < eps:
                if verbose:
                    print("  [k=%d] ||dx|| < eps=%.3f --> ostanov." % (k, eps))
                break
        else:
            if verbose:
                print("  [k=%d] phi_max = %.6f" % (k, max_phi))

        # Naibolee narushennoe ogranichenie
        i_viol = int(np.argmax(phis))
        g = grad_phi(x_k, i_viol + 1)
        phi_val = phis[i_viol]

        # Otsekayuschaya giperploskost s parametrom alpha (sm. dokazatelstvo)
        a_new = g / alpha
        b_new = (g @ x_k - phi_val) / alpha
        A_cuts.append(a_new)
        b_cuts.append(b_new)

        if verbose:
            print("       Narusheno phi_%d(x_k)=%.4f, dobavlena otchechka (alpha=%.2f)"
                  % (i_viol + 1, phi_val, alpha))

        x_prev = x_k.copy()
        res = solve_lp()
        if not res.success:
            if verbose:
                print("  [k=%d] Oshibka LP: %s" % (k, res.message))
            break
        x_k = res.x.copy()

        if verbose:
            print("       x_%d = %s,  f = %.5f"
                  % (k + 1, np.round(x_k, 4), c_obj @ x_k))

    f_opt = float(c_obj @ x_k)

    dual_vars = None
    try:
        dual_vars = res.ineqlin.marginals
    except AttributeError:
        pass

    return x_k, f_opt, dual_vars, np.array(A_cuts), np.array(b_cuts)


def print_dual(A_cuts, b_cuts, mu):
    """
    Dvoystvenная zadacha dlya finalnoy LP:
        Pryamaya:     min  c^T x         s.t. A x <= b
        Dvoystvenная: max  b^T mu        s.t. A^T mu = c,  mu <= 0

    Po teoreme silnoy dvoystvennosti (LP): c^T x* = b^T mu*.
    """
    print("\n  --- Dvoystvenная zadacha ---")
    print("  Pryamaya:     min c^T x    s.t. A_k x <= b_k")
    print("  Dvoystvenная: max b_k^T mu  s.t. A_k^T mu = c,  mu <= 0")
    if mu is not None:
        nz = np.where(np.abs(mu) > 1e-6)[0]
        print("  Nenulvye mu (aktivnye otchechki):")
        for idx in nz:
            print("    mu[%d] = %.6f" % (idx, mu[idx]))
        dual_obj = float(b_cuts @ mu)
        print("  Znachenie dvoystv. tseli b^T mu = %.5f  (= pryamaya tsel po silnoy dvoystvennosti)"
              % dual_obj)
    else:
        print("  (dvoystv. peremennye nedostupny v ustanovlennoy versii scipy)")


# =============================================================
#                         ГЛАВНАЯ ЧАСТЬ
# =============================================================
if __name__ == "__main__":
    sep = "=" * 65

    print(sep)
    print("  LAB. 5, CHAST 1 -- METOD OTSEKAYUSCHIKH GIPERPLOSKOSTEY")
    print(sep)
    print("  Zadacha:  min  x1 + 2*x2 + x3 + 2*x4")
    print("  phi_1:   x1^2 + x2^2 + x3^2 + x4^2 - 10 <= 0  [kvadratichnoe]")
    print("  phi_2:   x1^4 + x2^4 + x3^4 + x4^4 - 40 <= 0")
    print("  phi_3:   1.2*x1^4 + x2^4 + 1.1*x3^4 + x4^4 - 42 <= 0")
    print("  Tochnoe: x* = (-1, -2, -1, -2),  f* = -10")
    print()

    alpha = 1.0

    for eps in [0.1, 0.01, 0.001]:
        print(sep)
        print("  eps = %.3f,  alpha = %.1f" % (eps, alpha))
        print(sep)

        x_opt, f_opt, mu, A_fin, b_fin = cutting_plane(eps=eps, alpha=alpha, verbose=True)

        print()
        print("  REZULTAT:  x* ~= %s" % np.round(x_opt, 4))
        print("             f* ~= %.4f  (tochnoe: -10)" % f_opt)
        print("  Proverka ogranicheniy:")
        for i in range(1, 4):
            val = phi(x_opt, i)
            status = "<= 0 OK" if val <= 1e-3 else ">> NARUSHENO"
            print("    phi_%d(x*) = %.4f  %s" % (i, val, status))

        print_dual(A_fin, b_fin, mu)
        print()
