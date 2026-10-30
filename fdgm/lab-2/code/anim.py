"""
Визуализация счёта в реальном времени (задание 1б).

На экране одновременно показываются текущее численное решение v(x, t_n)
и установившееся решение w(x), к которому оно стремится.

Примеры запуска:

    python anim.py                      явная схема, M0, tau = tau0
    python anim.py --scheme implicit    неявная схема
    python anim.py --mult 1.02          явная схема за порогом устойчивости
    python anim.py --save evolution.gif сохранить анимацию в файл
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

import problem as pb
import schemes as sc


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--scheme', default='explicit',
                    choices=['explicit', 'implicit', 'symmetric'])
    ap.add_argument('--M', type=int, default=151,
                    help='число разбиений по x (по умолчанию найденное M0)')
    ap.add_argument('--mult', type=float, default=1.0,
                    help='множитель к шагу tau0 = h^2/(2 rho_max)')
    ap.add_argument('--fps', type=int, default=25)
    ap.add_argument('--frames', type=int, default=200,
                    help='сколько кадров показать за весь расчёт')
    ap.add_argument('--save', default=None, help='сохранить в GIF')
    return ap.parse_args()


def main():
    args = parse_args()
    T = pb.t_settle()
    M = args.M
    h = (pb.B - pb.A) / M
    tau = args.mult * sc.tau_stab(h)
    theta = {'explicit': 0.0, 'implicit': 1.0, 'symmetric': 0.5}[args.scheme]

    x = np.linspace(pb.A, pb.B, M + 1)
    dat = sc.Data(x)
    r = tau * dat.rho / h ** 2
    L0, LM = sc.boundary_rows(h)
    solver = None
    k0 = kM = 0.0
    if theta > 0:
        solver, (k0, kM) = sc.build_matrix(x, h, tau, theta)

    v = pb.phi(x)
    lap = np.zeros_like(v)
    nsteps = int(np.ceil(T / tau))
    every = max(1, nsteps // args.frames)

    state = {'v': v, 'n': 0}

    def step():
        v = state['v']
        t = state['n'] * tau
        t_new = t + tau
        if theta == 0.0:
            lap[1:-1] = r[1:-1] * (v[:-2] - 2 * v[1:-1] + v[2:])
            v = v + lap + tau * dat.f(t)
            v[0] = (dat.ga(t_new) - L0[1] * v[1] - L0[2] * v[2]) / L0[0]
            v[M] = (dat.gb(t_new) - LM[1] * v[M - 1] - LM[2] * v[M - 2]) / LM[0]
        else:
            rhs = v + tau * ((1.0 - theta) * dat.f(t) + theta * dat.f(t_new))
            if theta < 1.0:
                lap[1:-1] = r[1:-1] * (v[:-2] - 2 * v[1:-1] + v[2:])
                rhs += (1.0 - theta) * lap
            rhs[0] = dat.ga(t_new)
            rhs[M] = dat.gb(t_new)
            rhs[0] -= k0 * rhs[1]
            rhs[M] -= kM * rhs[M - 1]
            v = solver.solve(rhs)
        state['v'] = v
        state['n'] += 1

    xf = np.linspace(pb.A, pb.B, 2001)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(xf, pb.u_st(xf), 'k--', lw=1.6, label='установившееся $w(x)$')
    line, = ax.plot(x, v, lw=1.6, color='#c0392b', label='численное $v$')
    ax.set_xlabel('$x$')
    ax.set_ylabel('$u$')
    ax.grid(alpha=0.3)
    ax.legend(loc='upper right', fontsize=9)
    lim = 1.6 * np.max(np.abs(pb.phi(x)))
    ax.set_ylim(-lim, lim)
    title = ax.set_title('')

    def update(_):
        for _ in range(every):
            if state['n'] >= nsteps:
                break
            step()
        line.set_ydata(state['v'])
        title.set_text('%s, $M=%d$, $\\tau=%.3g$ ($%.2f\\,\\tau_0$), '
                       '$t=%.3f$ из %.3f'
                       % (args.scheme, M, tau, args.mult,
                          state['n'] * tau, T))
        return line, title

    anim = FuncAnimation(fig, update, frames=args.frames,
                         interval=1000 // args.fps, blit=False, repeat=False)
    if args.save:
        anim.save(args.save, writer='pillow', fps=args.fps)
        print('сохранено:', args.save)
    else:
        plt.show()


if __name__ == '__main__':
    main()
