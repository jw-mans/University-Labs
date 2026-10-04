"""
Расчеты по лабораторной: сходимость, ошибки округления, метод Рунге.

Печатает таблицы в консоль и строит графики в ../report/figures.
"""

import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import problem as pb
import schemes as sc

FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   '..', 'report', 'figures')
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({'font.size': 11, 'figure.dpi': 120,
                     'axes.grid': True, 'grid.alpha': 0.3})

# последовательность сеток: шаг уменьшается в 2 раза
NS = [20 * 2 ** j for j in range(14)]


def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight')
    plt.close(fig)


# отладка
def debug_quadratic():
    """
    На решении u* = x^2 схема O(h^2) должна давать машинную точность.
    """
    print('Отладка, u* = x^2')
    print('%8s %8s %12s' % ('n', 'схема', '||z_h||'))
    for n in (20, 40, 80):
        for order in (1, 2):
            x, v = sc.solve(n, order, f=pb.f_dbg,
                            ga=pb.gamma_a_dbg(), gb=pb.gamma_b_dbg())
            print('%8d %8s %12.3e'
                  % (n, 'O(h^%d)' % order, sc.err_c(v, pb.u_dbg(x))))


# графики решения
def fig_exact():
    x = np.linspace(pb.A, pb.B, 2001)
    fig, ax = plt.subplots(figsize=(7, 3.6))
    ax.plot(x, pb.u_exact(x), lw=1.6, color='#1f4e79')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u^*(x)$')
    savefig(fig, 'exact')


def fig_compare(n=40):
    xf = np.linspace(pb.A, pb.B, 2001)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    for ax, order in zip(axes, (1, 2)):
        x, v = sc.solve(n, order)
        ax.plot(xf, pb.u_exact(xf), lw=1.4, color='#1f4e79', label='$u^*$')
        ax.plot(x, v, 'o--', ms=3.5, lw=0.9, color='#c0392b', label='$v_h$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('$u$')
        ax.set_title('$n=%d$, $O(h^%d)$' % (n, order))
        ax.legend(fontsize=9)
    savefig(fig, 'compare')


# сходимость
def convergence(dtype=np.float64):
    """Погрешность обеих схем на всех сетках."""
    out = {1: [], 2: []}
    for n in NS:
        for order in (1, 2):
            x, v = sc.solve(n, order, dtype=dtype)
            out[order].append(sc.err_c(v, pb.u_exact(x)))
    return {k: np.array(val) for k, val in out.items()}


def orders(err):
    """Эмпирический порядок точности по двум соседним сеткам."""
    return [np.nan] + [np.log2(err[j - 1] / err[j])
                       for j in range(1, len(err))]


def fig_convergence(err):
    h = np.array([(pb.B - pb.A) / n for n in NS])
    mlh = -np.log10(h)
    fig, ax = plt.subplots(figsize=(7, 5))
    style = {1: ('o-', '#c0392b'), 2: ('s-', '#1f4e79')}
    for order in (1, 2):
        mlz = -np.log10(err[order])
        ax.plot(mlh, mlz, style[order][0], color=style[order][1], ms=4,
                label='схема $O(h^%d)$' % order)
        # эталонная прямая с наклоном k
        m = len(mlh) // 2
        c = mlz[m] - order * mlh[m]
        ax.plot(mlh, order * mlh + c, '--', lw=1.0, color='gray')
        ax.annotate('$k=%d$' % order, (mlh[-1], order * mlh[-1] + c),
                    textcoords='offset points', xytext=(-30, 6), color='gray')
    ax.set_xlabel(r'$-\lg h$')
    ax.set_ylabel(r'$-\lg \|z_h\|$')
    ax.legend()
    savefig(fig, 'convergence')


# ошибки округления
def fig_roundoff(err64, err32):
    h = np.array([(pb.B - pb.A) / n for n in NS])
    mlh = -np.log10(h)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for ax, order in zip(axes, (1, 2)):
        ax.plot(mlh, -np.log10(err64[order]), 'o-', ms=4, color='#1f4e79',
                label='double (float64)')
        ax.plot(mlh, -np.log10(err32[order]), 's-', ms=4, color='#c0392b',
                label='single (float32)')
        ax.set_xlabel(r'$-\lg h$')
        ax.set_title('схема $O(h^%d)$' % order)
        ax.legend(fontsize=9)
    axes[0].set_ylabel(r'$-\lg \|z_h\|$')
    savefig(fig, 'roundoff')


def min_nodes(order, eps, n_max=200000):
    """
    Наименьшее n, при котором ||z_h|| <= eps.

    Погрешность монотонно убывает с ростом n, поэтому нужное n ищется
    двоичным поиском, а не выбирается из последовательности удвоений.
    """
    def err(n):
        x, v = sc.solve(n, order)
        return sc.err_c(v, pb.u_exact(x))

    if err(n_max) > eps:
        return None, None

    lo, hi = 2, n_max # err(lo) > eps >= err(hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if err(mid) <= eps:
            hi = mid
        else:
            lo = mid
    return hi, err(hi)


# метод Рунге
def runge(order, ns):
    """
    Оценка погрешности по Рунге и уточнение решения.

    Возвращает строки (n, ||R||, ||z_{h/2}||, ||u - v_уточн||).
    """
    k = order
    sols = {n: sc.solve(n, order) for n in ns}
    rows = []
    for n in ns[:-1]:
        x1, v1 = sols[n]
        v2 = sols[2 * n][1][::2] # решение на мелкой сетке в общих узлах
        R = (v2 - v1) / (2 ** k - 1) # оценка погрешности v_{h/2}
        ue = pb.u_exact(x1)
        rows.append((n,
                     np.max(np.abs(R)),
                     np.max(np.abs(ue - v2)),
                     np.max(np.abs(ue - (v2 + R)))))
    return rows


def even_check(order, ns, node='a'):
    """
    Проверка, какие степени h входят в разложение погрешности.

    Если z(h) = c_k h^k + c_{k+1} h^{k+1} + ..., то комбинация
    d(h) = z(h) - 2^k z(h/2) убирает главный член. Тогда отношение
    d(h)/d(h/2) равно 2^{k+1}, если следующий член нечётный, и
    2^{k+2}, если он чётный (нечётного нет).

    Погрешность берётся в одном узле, чтобы не переключаться между
    точками максимума. Возвращает строки (n, d, d/d2).
    """
    z = {}
    for n in ns:
        x, v = sc.solve(n, order)
        i = {'a': 0, 'mid': len(x) // 2, 'b': len(x) - 1}[node]
        z[n] = pb.u_exact(x[i]) - v[i]
    d = {n: z[n] - 2 ** order * z[2 * n] for n in ns[:-1]}
    return [(n, d[n], d[n] / d[2 * n]) for n in ns[:-2]]


def fig_runge(rows1, rows2):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, rows, order, kref in ((axes[0], rows1, 1, 2),
                                  (axes[1], rows2, 2, 4)):
        h = (pb.B - pb.A) / np.array([r[0] for r in rows])
        mlh = -np.log10(h)
        z2 = np.array([r[2] for r in rows])
        zr = np.array([r[3] for r in rows])
        ax.plot(mlh, -np.log10(z2), 'o-', ms=4, color='#1f4e79',
                label=r'$\|z_{h/2}\|$')
        ax.plot(mlh, -np.log10(zr), 's-', ms=4, color='#c0392b',
                label='после уточнения')
        m = len(mlh) // 2
        for kk, base in ((order, z2), (kref, zr)):
            c = -np.log10(base[m]) - kk * mlh[m]
            ax.plot(mlh, kk * mlh + c, '--', lw=1.0, color='gray')
            ax.annotate('$k=%d$' % kk, (mlh[-1], kk * mlh[-1] + c),
                        textcoords='offset points', xytext=(-30, 6),
                        color='gray', fontsize=9)
        ax.set_xlabel(r'$-\lg h$')
        ax.set_ylabel(r'$-\lg \|z\|$')
        ax.set_title('схема $O(h^%d)$' % order)
        ax.legend(fontsize=9)
    savefig(fig, 'runge')


# вывод
def main():
    debug_quadratic()
    fig_exact()
    fig_compare()

    # задания 2-4: сходимость
    e64 = convergence(np.float64)
    fig_convergence(e64)
    k1, k2 = orders(e64[1]), orders(e64[2])

    print('\nСходимость (двойная точность)')
    print('%8s %10s %12s %6s %12s %6s'
          % ('n', 'h', '||z|| O(h)', 'k', '||z|| O(h^2)', 'k'))
    for j, n in enumerate(NS):
        print('%8d %10.2e %12.2e %6.2f %12.2e %6.2f'
              % (n, (pb.B - pb.A) / n, e64[1][j], k1[j], e64[2][j], k2[j]))

    # задание 5: ошибки округления
    e32 = convergence(np.float32)
    fig_roundoff(e64, e32)

    print('\nОдинарная и двойная точность')
    print('%8s %12s %12s %12s %12s'
          % ('n', 'O(h) dbl', 'O(h) sgl', 'O(h^2) dbl', 'O(h^2) sgl'))
    for j, n in enumerate(NS):
        print('%8d %12.2e %12.2e %12.2e %12.2e'
              % (n, e64[1][j], e32[1][j], e64[2][j], e32[2][j]))

    print('\nМинимум погрешности в одинарной точности')
    for order in (1, 2):
        j = int(np.argmin(e32[order]))
        print('  O(h^%d): n = %d, ||z|| = %.2e' % (order, NS[j], e32[order][j]))

    # задание 6: число узлов для заданной точности
    print('\nЧисло разбиений для заданной точности (двойная точность)')
    print('%10s %10s %12s %10s %12s'
          % ('eps', 'n, O(h)', '||z||', 'n, O(h^2)', '||z||'))
    for eps in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
        n1, z1 = min_nodes(1, eps)
        n2, z2 = min_nodes(2, eps)
        print('%10.0e %10s %12s %10s %12s'
              % (eps,
                 n1 or '-', '%.2e' % z1 if n1 else '-',
                 n2 or '-', '%.2e' % z2 if n2 else '-'))

    # допзадание: метод Рунге
    r1 = runge(1, [20 * 2 ** j for j in range(10)])
    r2 = runge(2, [20 * 2 ** j for j in range(8)])
    fig_runge(r1, r2)

    for order, rows in ((1, r1), (2, r2)):
        print('\nМетод Рунге, схема O(h^%d)' % order)
        print('%8s %12s %12s %8s %12s %7s'
              % ('n', '||R||', '||z_h/2||', 'R/z', 'после', 'порядок'))
        prev = None
        for n, R, z2, zr in rows:
            k = '-' if prev is None else '%.2f' % np.log2(prev / zr)
            print('%8d %12.2e %12.2e %8.3f %12.2e %7s'
                  % (n, R, z2, R / z2, zr, k))
            prev = zr

    print('\nКакие степени h входят в погрешность (узел x = a)')
    print('%8s %8s %14s %8s' % ('схема', 'n', 'd(h)', 'd/d2'))
    for order in (1, 2):
        for n, d, r in even_check(order, [40 * 2 ** j for j in range(6)]):
            print('%8s %8d %14.3e %8.2f' % ('O(h^%d)' % order, n, d, r))

    print('\nГрафики:', os.path.normpath(FIG))


if __name__ == '__main__':
    main()
