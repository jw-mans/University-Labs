"""
ЛР 2. Постановка задачи.

Краевая задача третьего рода для параболического уравнения:

    du/dt - rho(x) d^2u/dx^2 = f(x, t),   a <= x <= b,  0 <= t <= T,
    0 < rho_min <= rho(x) <= rho_max,     rho(x) != const,
    u(x, 0) = phi(x),

краевые условия ставятся так же, как в работе 1 (третьего рода):

    -alpha_a du/dx(a, t) + beta_a u(a, t) = gamma_a(t),
     alpha_b du/dx(b, t) + beta_b u(b, t) = gamma_b(t).

Точное решение строится в виде u* = v(x, t) w(x), где v -> 1 при t -> oo:

    v(x, t) = p(x) exp(-k t) + 1,   k > 0,
    w(x)    -- решение из работы 1 (установившееся состояние).

Тогда u*(x, t) -> w(x) при t -> oo, а w удовлетворяет стационарной задаче
rho(x) w'' = g(x) с теми же краевыми условиями. Правые части f, phi,
gamma_a, gamma_b вычисляются подстановкой u* в уравнение.
"""

import numpy as np

# отрезок по x
A, B = 0.0, 2.0

# граничные коэффициенты (те же, что в работе 1)
ALPHA_A, BETA_A = 1.0, 1.0
ALPHA_B, BETA_B = 1.0, 2.0

# скорость выхода на установившийся режим
K = 3.0

# порог установления: ||u(.,t) - u_st|| / ||u_st|| <= DELTA
DELTA = 1e-3

RHO_MIN, RHO_MAX = 1.0, 2.0


def rho(x):
    """Коэффициент теплопроводности. Не константа, 1 <= rho <= 2.

    Максимум rho = 2 достигается при x = 1, то есть во внутренней точке
    отрезка. Это важно: граница устойчивости явной схемы определяется
    максимумом rho по тем узлам, которые пересчитываются по явной
    формуле, то есть по внутренним. При максимуме на краю оценка
    tau <= h^2/(2 rho_max) была бы завышенной.
    """
    return 1.0 + np.sin(np.pi * x / 2) ** 2


# установившееся решение w(x) (из работы 1)
def w(x):
    return x * np.sin(2 * np.pi * x) + np.exp(-x)


def dw(x):
    return (np.sin(2 * np.pi * x)
            + 2 * np.pi * x * np.cos(2 * np.pi * x)
            - np.exp(-x))


def d2w(x):
    return (4 * np.pi * np.cos(2 * np.pi * x)
            - 4 * np.pi ** 2 * x * np.sin(2 * np.pi * x)
            + np.exp(-x))


# множитель p(x) в v(x, t)
def p(x):
    return 1.0 + x / 2


# p' = 1/2, p'' = 0, поэтому
def pw(x):
    return p(x) * w(x)


def d_pw(x):
    return 0.5 * w(x) + p(x) * dw(x)


def d2_pw(x):
    return dw(x) + p(x) * d2w(x)


# точное решение u*
def u_exact(x, t):
    return pw(x) * np.exp(-K * t) + w(x)


def du_dx(x, t):
    return d_pw(x) * np.exp(-K * t) + dw(x)


def d2u_dx2(x, t):
    return d2_pw(x) * np.exp(-K * t) + d2w(x)


def du_dt(x, t):
    return -K * pw(x) * np.exp(-K * t)


# согласованные данные задачи
def f_star(x, t):
    """f = du/dt - rho d^2u/dx^2."""
    return du_dt(x, t) - rho(x) * d2u_dx2(x, t)


def phi(x):
    """Начальное условие u*(x, 0)."""
    return u_exact(x, 0.0)


def gamma_a(t):
    return -ALPHA_A * du_dx(A, t) + BETA_A * u_exact(A, t)


def gamma_b(t):
    return ALPHA_B * du_dx(B, t) + BETA_B * u_exact(B, t)


def u_st(x):
    """Установившееся решение, к которому стремится u* при t -> oo."""
    return w(x)


# момент установления
def t_settle(delta=DELTA, n=200001):
    """Момент ust, начиная с которого ||u(.,t)-u_st||/||u_st|| <= delta.

    Так как u(x,t) - u_st(x) = p(x) w(x) exp(-k t), неравенство решается
    точно:  t >= (1/k) ln( ||p w|| / (delta ||w||) ).
    Нормы несеточные, считаются по мелкой сетке.
    """
    x = np.linspace(A, B, n)
    return np.log(np.max(np.abs(pw(x)))
                  / (delta * np.max(np.abs(w(x))))) / K
