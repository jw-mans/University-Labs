"""
ЛР 1. Постановка задачи.

Краевая задача третьего рода:

    (A u)(x) := -u'' + p(x) u' + q(x) u = f(x),      x in [a, b]
    (A u)(a) := -alpha_a u'(a) + beta_a u(a) = gamma_a
    (A u)(b) :=  alpha_b u'(b) + beta_b u(b) = gamma_b

Задаем точное решение u*, вычисляем по нему правые части
f* = (A u*)(x), gamma*_a = (A u*)(a), gamma*_b = (A u*)(b).
Тогда точное решение полученной задачи -- само u*, и погрешность
численного решения можно посчитать в каждом узле.
"""

import numpy as np

# отрезок
A, B = 0.0, 2.0

# граничные коэффициенты
ALPHA_A, BETA_A = 1.0, 1.0     # -alpha_a u'(a) + beta_a u(a) = gamma_a
ALPHA_B, BETA_B = 1.0, 2.0     #  alpha_b u'(b) + beta_b u(b) = gamma_b


def p(x):
    """Коэффициент при u'. Не константа."""
    return 2.0 + np.cos(np.pi * x)


def q(x):
    """Коэффициент при u. Не константа, q >= 1 > 0."""
    return 1.0 + x ** 2


# точное решение u* = x sin(2 pi x) + exp(-x)
def u_exact(x):
    return x * np.sin(2 * np.pi * x) + np.exp(-x)


def du_exact(x):
    return (np.sin(2 * np.pi * x)
            + 2 * np.pi * x * np.cos(2 * np.pi * x)
            - np.exp(-x))


def d2u_exact(x):
    return (4 * np.pi * np.cos(2 * np.pi * x)
            - 4 * np.pi ** 2 * x * np.sin(2 * np.pi * x)
            + np.exp(-x))


# согласованные правые части
def f_star(x):
    return -d2u_exact(x) + p(x) * du_exact(x) + q(x) * u_exact(x)


def gamma_a_star():
    return -ALPHA_A * du_exact(A) + BETA_A * u_exact(A)


def gamma_b_star():
    return ALPHA_B * du_exact(B) + BETA_B * u_exact(B)


# отладочное решение u* = x^2
# схема O(h^2) должна воспроизводить его точно
def u_dbg(x):
    return x ** 2


def du_dbg(x):
    return 2 * x


def d2u_dbg(x):
    return np.full_like(np.asarray(x, dtype=float), 2.0)


def f_dbg(x):
    return -d2u_dbg(x) + p(x) * du_dbg(x) + q(x) * u_dbg(x)


def gamma_a_dbg():
    return -ALPHA_A * du_dbg(A) + BETA_A * u_dbg(A)


def gamma_b_dbg():
    return ALPHA_B * du_dbg(B) + BETA_B * u_dbg(B)
