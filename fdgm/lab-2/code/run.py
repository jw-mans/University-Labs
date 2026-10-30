"""
Расчёты по лабораторной 2: явная, неявная и симметричная схемы.

Печатает таблицы в консоль и строит графики в ../report/figures.
"""

import os
import time

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

EPS0 = 5e-3                     # заданная точность
DTAU_REL = 0.02                 # относительная добавка к tau0 в задании 1в
NCTRL = 101                     # число контрольных моментов времени
COL = {'explicit': '#c0392b', 'implicit': '#1f4e79', 'symmetric': '#1e8449'}
RUS = {'explicit': 'явная', 'implicit': 'неявная',
       'symmetric': 'симметричная'}


def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight')
    plt.close(fig)


def eps_explicit(M, T):
    """Погрешность явной схемы при максимально допустимом шаге по времени."""
    h = (pb.B - pb.A) / M
    return sc.solve_pde(M, sc.tau_stab(h), 'explicit', T)['eps']


# Задание 1а: число разбиений M0
def find_M0(T, eps0=EPS0, m_max=800):
    """
    Наименьшее M, при котором eps(M) < eps0.

    Погрешность убывает с ростом M, поэтому сначала сетка удваивается до
    выполнения неравенства (это даёт вилку), а затем M0 уточняется
    двоичным поиском. Удвоение снизу, а не проверка сразу при большом M,
    экономит счёт: стоимость одного прогона растёт как M^3 (число шагов
    по времени само пропорционально M^2).
    """
    lo, hi = 8, 16
    while eps_explicit(hi, T) >= eps0:
        lo, hi = hi, 2 * hi
        if hi > m_max:
            raise RuntimeError('не хватает m_max')
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if eps_explicit(mid, T) < eps0:
            hi = mid
        else:
            lo = mid
    return hi, eps_explicit(hi, T), eps_explicit(hi - 1, T)


# Задание 4 и 5: максимальный шаг по времени при заданной точности
def find_tau_max(M, scheme, T, tau_lo, eps0=EPS0, grow=2.0, key='eps'):
    """
    Максимальное tau, при котором ещё eps(tau) <= eps0.

    Сначала шаг увеличивается в grow раз до нарушения неравенства,
    затем граница уточняется двоичным поиском. key выбирает, какая из
    двух оценок погрешности используется: 'eps' -- по сеточной области
    Omega_h,tau (как в задании), 'eps_ctrl' -- по фиксированному набору
    контрольных моментов времени.
    """
    ctrl = np.linspace(0.0, T, NCTRL) if key == 'eps_ctrl' else None

    def eps(tau):
        return sc.solve_pde(M, tau, scheme, T, ctrl=ctrl)[key]

    lo = tau_lo
    hi = tau_lo * grow
    while eps(hi) <= eps0:
        lo, hi = hi, hi * grow
        if hi > T:
            return T, eps(T)
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if eps(mid) <= eps0:
            lo = mid
        else:
            hi = mid
        if (hi - lo) / lo < 1e-3:
            break
    return lo, eps(lo)


# Графики
def fig_exact(T):
    x = np.linspace(pb.A, pb.B, 2001)
    times = [0.0, 0.1, 0.25, 0.5, 1.0, T]
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, len(times)))
    for c, t in zip(cmap, times):
        ax.plot(x, pb.u_exact(x, t), lw=1.3, color=c, label='$t=%.2f$' % t)
    ax.plot(x, pb.u_st(x), 'k--', lw=1.6, label='$w(x)$')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u^*(x,t)$')
    ax.legend(fontsize=8, ncol=2)
    savefig(fig, 'exact')


def fig_rho():
    x = np.linspace(pb.A, pb.B, 1001)
    fig, ax = plt.subplots(figsize=(6, 2.6))
    ax.plot(x, pb.rho(x), lw=1.6, color='#1f4e79')
    ax.set_xlabel('$x$')
    ax.set_ylabel(r'$\rho(x)$')
    ax.set_ylim(0.9, 2.1)
    savefig(fig, 'rho')


def fig_convergence(T, Ms, errs, M0):
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    h = np.array([(pb.B - pb.A) / m for m in Ms])
    ax.plot(-np.log10(h), -np.log10(errs), 'o-', ms=4, color=COL['explicit'],
            label='явная схема')
    m = len(h) // 2
    c = -np.log10(errs[m]) - 2 * (-np.log10(h[m]))
    ax.plot(-np.log10(h), 2 * (-np.log10(h)) + c, '--', lw=1.0,
            color='gray', label='эталонная прямая, $k=2$')
    ax.axhline(-np.log10(EPS0), color='#8e44ad', lw=1.0, ls=':')
    ax.annotate(r'$\varepsilon_0$', (-np.log10(h[0]), -np.log10(EPS0)),
                textcoords='offset points', xytext=(4, 4), color='#8e44ad')
    ax.set_xlabel(r'$-\lg h$')
    ax.set_ylabel(r'$-\lg \varepsilon$')
    ax.legend(fontsize=9)
    savefig(fig, 'convergence')


def fig_frames(res, name, title):
    """Кадры численного решения и «висящее» установившееся решение."""
    xf = np.linspace(pb.A, pb.B, 2001)
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    cmap = plt.cm.viridis(np.linspace(0.05, 0.85, len(res['frames'])))
    for c, (t, v) in zip(cmap, res['frames']):
        ax.plot(res['x'], v, lw=1.2, color=c, label='$t=%.2f$' % t)
    ax.plot(xf, pb.u_st(xf), 'k--', lw=1.6, label='$w(x)$')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$v$')
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, ncol=2)
    savefig(fig, name)


def fig_instability(res_bad, res_ok, tau0, dtau, M0):
    """Слева -- рост погрешности во времени, справа -- пилообразный
    профиль погрешности перед разрушением счёта."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    t_end = res_bad['t_end']

    ax = axes[0]
    ax.semilogy(res_ok['thist'], res_ok['ehist'], lw=1.4,
                color=COL['explicit'], label=r'$\tau=\tau_0$')
    ax.semilogy(res_bad['thist'], res_bad['ehist'], lw=1.4,
                color='#8e44ad', label=r'$\tau=\tau_0+\Delta\tau$')
    # теоретический рост: спектральный радиус матрицы перехода.
    # Прямая проводится через последнюю точку назад по времени.
    lam = sc.amplification(M0, tau0 + dtau)
    e, th = res_bad['ehist'], res_bad['thist']
    tt = np.linspace(0.45 * th[-1], th[-1], 50)
    ax.semilogy(tt, e[-1] * lam ** ((tt - th[-1]) / (tau0 + dtau)),
                '--', lw=1.2, color='gray',
                label=r'$\lambda_{\max}^{\,n}$, $\lambda_{\max}=%.5f$' % lam)
    ax.set_xlim(0, 1.15 * t_end)
    ax.set_ylim(1e-3, 1e3)
    ax.set_xlabel('$t$')
    ax.set_ylabel(r'$\max_i|u(x_i,t_n)-v_i^n|$')
    ax.set_title('рост погрешности', fontsize=10)
    ax.legend(fontsize=8)

    ax = axes[1]
    x = res_bad['x']
    z = res_bad['v'] - pb.u_exact(x, t_end)
    j = int(np.argmax(np.abs(z)))
    lo, hi = max(0, j - 12), min(x.size, j + 13)
    ax.plot(x[lo:hi], z[lo:hi], 'o-', ms=3.5, lw=1.0, color='#8e44ad')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$v_i^n-u(x_i,t_n)$')
    ax.set_title(r'профиль погрешности при $t=%.4f$ (увеличено)' % t_end,
                 fontsize=10)
    savefig(fig, 'instability')


def fig_eps_tau(curves, tau0, eps0=EPS0):
    """Зависимость погрешности от шага по времени для обеих схем."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, scheme in zip(axes, ('implicit', 'symmetric')):
        taus, eps, ectrl, tstar, tctrl = curves[scheme]
        rel = np.array(taus) / tau0
        ax.plot(rel, eps, 'o-', ms=4, color=COL[scheme],
                label=r'$\varepsilon$ по $\Omega_{h\tau}$')
        ax.plot(rel, ectrl, 's--', ms=4, color='#8e44ad',
                label=r'$\varepsilon$ по контрольным моментам')
        ax.axhline(eps0, color='k', ls='--', lw=1.0)
        ax.axvline(tstar / tau0, color=COL[scheme], ls=':', lw=1.2)
        ax.axvline(tctrl / tau0, color='#8e44ad', ls=':', lw=1.2)
        if scheme == 'symmetric':
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_ylim(4e-3, 5e-2)
        ax.set_xlabel(r'$\tau/\tau_0$')
        ax.set_ylabel(r'$\varepsilon$')
        ax.set_title(RUS[scheme] + ' схема', fontsize=10)
        ax.legend(fontsize=8)
    savefig(fig, 'eps_tau')


def main():
    T = pb.t_settle()
    x = np.linspace(pb.A, pb.B, 2001)

    print('Параметры задачи')
    print('  отрезок [%g, %g], k = %g, delta = %g' % (pb.A, pb.B, pb.K,
                                                      pb.DELTA))
    print('  rho: от %.3f до %.3f' % (pb.rho(x).min(), pb.rho(x).max()))
    print('  момент установления ust = T = %.4f' % T)
    print('  заданная точность eps0 = %g' % EPS0)

    fig_rho()
    fig_exact(T)

    # задание 1а 
    print('\nЗадание 1а. Поиск числа разбиений')
    t_begin = time.perf_counter()
    M0, eps_M0, eps_M0m1 = find_M0(T)
    h0 = (pb.B - pb.A) / M0
    tau0 = sc.tau_stab(h0)
    print('  M0 = %d,  h0 = %.6f,  tau0 = %.6e' % (M0, h0, tau0))
    print('  eps(M0)   = %.4e  <  eps0 = %.1e' % (eps_M0, EPS0))
    print('  eps(M0-1) = %.4e  >  eps0' % eps_M0m1)
    print('  (поиск занял %.1f с)' % (time.perf_counter() - t_begin))

    Ms = [20, 30, 40, 60, 80, 120, 160, 240]
    errs = []
    print('\n%8s %10s %12s %10s %12s' % ('M', 'h', 'tau_0', 'N', 'eps'))
    for m in Ms:
        h = (pb.B - pb.A) / m
        res = sc.solve_pde(m, sc.tau_stab(h), 'explicit', T)
        errs.append(res['eps'])
        print('%8d %10.5f %12.3e %10d %12.4e'
              % (m, h, sc.tau_stab(h), res['nsteps'], res['eps']))
    fig_convergence(T, Ms, np.array(errs), M0)

    snaps = [0.05, 0.15, 0.3, 0.6, 1.2, T]

    # задание 1б
    print('\nЗадание 1б. Счёт явной схемой при M0, tau0')
    res_exp = sc.solve_pde(M0, tau0, 'explicit', T, snapshots=snaps,
                           norm_every=50)
    print('  N = %d шагов, eps = %.4e' % (res_exp['nsteps'], res_exp['eps']))
    fig_frames(res_exp, 'explicit',
               r'явная схема, $M_0=%d$, $\tau=\tau_0$' % M0)

    # задание 1в
    dtau = DTAU_REL * tau0
    print('\nЗадание 1в. Явная схема при tau = tau0 + dtau, dtau = %.0f%% tau0'
          % (100 * DTAU_REL))
    res_bad = sc.solve_pde(M0, tau0 + dtau, 'explicit', T,
                           norm_every=20, blow_at=1e2)
    print('  счёт разрушился: %s' % ('да' if res_bad['blown'] else 'нет'))
    print('  дошёл до t = %.4f из %.4f' % (res_bad['t_end'], T))
    fig_instability(res_bad, res_exp, tau0, dtau, M0)

    print('\n  порог устойчивости: спектральный радиус матрицы перехода'
          ' и результат счёта')
    print('  %10s %16s %14s' % ('tau/tau0', 'спектр. радиус', 'счёт'))
    for mult in (0.95, 0.99, 1.00, 1.01, 1.02, 1.05):
        rr = sc.solve_pde(M0, mult * tau0, 'explicit', T, track_error=False)
        print('  %10.2f %16.6f %14s'
              % (mult, sc.amplification(M0, mult * tau0),
                 'разрушение' if rr['blown'] else 'устойчиво'))
    lo, hi = 0.9, 1.2
    for _ in range(40):
        mid = 0.5 * (lo + hi)
        if sc.amplification(M0, mid * tau0) <= 1.0:
            lo = mid
        else:
            hi = mid
    print('  порог по спектральному радиусу: tau/tau0 = %.5f' % lo)

    # задание 2
    print('\nЗадание 2. Неявная схема')
    res_imp = sc.solve_pde(M0, tau0, 'implicit', T, snapshots=snaps)
    print('  tau = tau0:        eps = %.4e, N = %d'
          % (res_imp['eps'], res_imp['nsteps']))
    res_imp_bad = sc.solve_pde(M0, tau0 + dtau, 'implicit', T,
                               snapshots=snaps)
    print('  tau = tau0 + dtau: eps = %.4e, N = %d, разрушение: %s'
          % (res_imp_bad['eps'], res_imp_bad['nsteps'],
             'да' if res_imp_bad['blown'] else 'нет'))
    fig_frames(res_imp_bad, 'implicit',
               r'неявная схема, $M_0=%d$, $\tau=\tau_0+\Delta\tau$' % M0)

    # задание 4 и 5.3
    ctrl = np.linspace(0.0, T, NCTRL)
    curves = {}
    summary = {}
    for scheme in ('implicit', 'symmetric'):
        print('\nМаксимальный шаг по времени, схема: %s' % RUS[scheme])
        t_begin = time.perf_counter()
        tstar, eps_star = find_tau_max(M0, scheme, T, tau0)
        tctrl, eps_c = find_tau_max(M0, scheme, T, tau0, key='eps_ctrl')
        print('  по Omega_h,tau:        tau* = %.4e = %.1f tau0, eps = %.4e'
              % (tstar, tstar / tau0, eps_star))
        print('  по контрольным точкам: tau* = %.4e = %.1f tau0, eps = %.4e'
              % (tctrl, tctrl / tau0, eps_c))
        print('  (поиск занял %.1f с)' % (time.perf_counter() - t_begin))

        taus = np.unique(np.concatenate([
            np.geomspace(tau0, tstar, 12), [tstar, tctrl]]))
        res = [sc.solve_pde(M0, t, scheme, T, ctrl=ctrl) for t in taus]
        eps = [r['eps'] for r in res]
        ectrl = [r['eps_ctrl'] for r in res]
        curves[scheme] = (taus, eps, ectrl, tstar, tctrl)
        summary[scheme] = (tstar, eps_star, tctrl, eps_c)
        print('  %14s %9s %7s %12s %12s'
              % ('tau', 'tau/tau0', 'N', 'eps', 'eps_контр'))
        for t, r in zip(taus, res):
            print('  %14.4e %9.2f %7d %12.4e %12.4e'
                  % (t, t / tau0, r['nsteps'], r['eps'], r['eps_ctrl']))

    fig_eps_tau(curves, tau0)

    res_imp_star = sc.solve_pde(M0, summary['implicit'][0], 'implicit', T,
                                snapshots=snaps)
    fig_frames(res_imp_star, 'implicit_star',
               r'неявная схема, $M_0=%d$, $\tau=\tau^*$' % M0)
    res_sym_star = sc.solve_pde(M0, summary['symmetric'][2], 'symmetric', T,
                                snapshots=snaps)
    fig_frames(res_sym_star, 'symmetric_star',
               r'симметричная схема, $M_0=%d$, $\tau=\tau^{**}$' % M0)

    # задание 5.4: трудоёмкость
    print('\nЗадание 5.4. Трудоёмкость при M0 = %d' % M0)
    plan = [('explicit', tau0), ('implicit', summary['implicit'][2]),
            ('symmetric', summary['symmetric'][2])]
    print('%14s %12s %8s %9s %10s %12s'
          % ('схема', 'tau', 'tau/tau0', 'N', 'время, с', 'eps'))
    for scheme, tau in plan:
        t_begin = time.perf_counter()
        r = sc.solve_pde(M0, tau, scheme, T, track_error=False)
        dt = time.perf_counter() - t_begin
        e = sc.solve_pde(M0, tau, scheme, T)['eps']
        print('%14s %12.4e %8.1f %9d %10.2f %12.4e'
              % (RUS[scheme], tau, tau / tau0, r['nsteps'], dt, e))

    print('\nГрафики:', os.path.normpath(FIG))


if __name__ == '__main__':
    main()
