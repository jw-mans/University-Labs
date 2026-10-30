"""
Разностные схемы для параболического уравнения: явная, неявная,
симметричная (Кранка -- Николсон). Прогонка с заранее вычисленными
прогоночными коэффициентами.
"""

import numpy as np

import problem as pb


# Прогонка
class Thomas:
    """
    Прогонка для трёхдиагональной матрицы с постоянными во времени
    коэффициентами: прямой ход факторизации выполняется один раз.

    lower[i] -- коэффициент при y[i-1] в строке i,
    upper[i] -- коэффициент при y[i+1] в строке i.
    """

    def __init__(self, lower, diag, upper):
        n = diag.size
        self.n = n
        self.m = np.zeros(n)          # множители прямого хода
        self.c = upper.copy()
        self.b = diag.copy()
        self.a = lower.copy()
        for i in range(1, n):
            self.m[i] = self.a[i] / self.b[i - 1]
            self.b[i] -= self.m[i] * self.c[i - 1]

    def solve(self, rhs):
        d = rhs.copy()
        for i in range(1, self.n):
            d[i] -= self.m[i] * d[i - 1]
        y = np.empty(self.n)
        y[-1] = d[-1] / self.b[-1]
        for i in range(self.n - 2, -1, -1):
            y[i] = (d[i] - self.c[i] * y[i + 1]) / self.b[i]
        return y


# Аппроксимация краевых условий (трёхточечная, второй порядок по h)
def boundary_rows(h):
    """
    Коэффициенты краевых уравнений.

    du/dx(a) = (-3 v_0 + 4 v_1 - v_2) / (2h) + O(h^2),
    du/dx(b) = ( 3 v_M - 4 v_{M-1} + v_{M-2}) / (2h) + O(h^2).

    Подстановка в краевые условия даёт

        L0 = (c00, c01, c02) для (v_0, v_1, v_2),
        LM = (cM0, cM1, cM2) для (v_M, v_{M-1}, v_{M-2}).
    """
    aa, ba = pb.ALPHA_A, pb.BETA_A
    ab, bb = pb.ALPHA_B, pb.BETA_B
    L0 = (1.5 * aa / h + ba, -2.0 * aa / h, 0.5 * aa / h)
    LM = (1.5 * ab / h + bb, -2.0 * ab / h, 0.5 * ab / h)
    return L0, LM


# Сборка матрицы неявной и симметричной схем
def build_matrix(x, h, tau, theta):
    """
    Матрица схемы с весом theta.

    theta = 0   -- явная,
    theta = 1   -- неявная,
    theta = 1/2 -- симметричная.

    Возвращает объект Thomas для системы относительно v^{n+1}.
    Краевые строки трёхточечные, поэтому лишние элементы (v_2 в строке 0
    и v_{M-2} в строке M) исключаются вычитанием соседней строки.
    """
    M = x.size - 1
    r = tau * pb.rho(x) / h ** 2

    lo = np.zeros(M + 1)
    di = np.ones(M + 1)
    up = np.zeros(M + 1)

    i = np.arange(1, M)
    lo[i] = -theta * r[i]
    di[i] = 1.0 + 2.0 * theta * r[i]
    up[i] = -theta * r[i]

    L0, LM = boundary_rows(h)
    # строка 0: c00 v_0 + c01 v_1 + c02 v_2 = ...
    # исключаем v_2, вычитая строку 1, умноженную на c02 / up[1]
    k0 = L0[2] / up[1]
    di[0] = L0[0] - k0 * lo[1]
    up[0] = L0[1] - k0 * di[1]
    # строка M: cM0 v_M + cM1 v_{M-1} + cM2 v_{M-2} = ...
    kM = LM[2] / lo[M - 1]
    di[M] = LM[0] - kM * up[M - 1]
    lo[M] = LM[1] - kM * di[M - 1]

    return Thomas(lo, di, up), (k0, kM)


# Шаг по времени
class Data:
    """
    Предвычисленные сеточные массивы.

    Все данные задачи зависят от t только через множитель exp(-k t):

        u*(x, t) = pw(x) e^{-kt} + w(x),
        f (x, t) = F1(x) e^{-kt} + F0(x),
        gamma(t) = g1 e^{-kt} + g0,

    поэтому функции от x считаются один раз, а на каждом шаге остаётся
    умножение на скаляр. Без этого основное время уходит на пересчёт
    синусов и экспонент.
    """

    def __init__(self, x):
        self.x = x
        self.rho = pb.rho(x)
        self.W = pb.w(x)
        self.PW = pb.pw(x)
        self.F1 = -pb.K * self.PW - self.rho * pb.d2_pw(x)
        self.F0 = -self.rho * pb.d2w(x)
        self.ga1 = -pb.ALPHA_A * pb.d_pw(pb.A) + pb.BETA_A * pb.pw(pb.A)
        self.ga0 = -pb.ALPHA_A * pb.dw(pb.A) + pb.BETA_A * pb.w(pb.A)
        self.gb1 = pb.ALPHA_B * pb.d_pw(pb.B) + pb.BETA_B * pb.pw(pb.B)
        self.gb0 = pb.ALPHA_B * pb.dw(pb.B) + pb.BETA_B * pb.w(pb.B)

    def f(self, t):
        return self.F1 * np.exp(-pb.K * t) + self.F0

    def u(self, t):
        return self.PW * np.exp(-pb.K * t) + self.W

    def ga(self, t):
        return self.ga1 * np.exp(-pb.K * t) + self.ga0

    def gb(self, t):
        return self.gb1 * np.exp(-pb.K * t) + self.gb0


def solve_pde(M, tau, scheme, T, track_error=True, snapshots=None,
              norm_every=0, ctrl=None, blow_at=1e8):
    """
    Проинтегрировать задачу до момента T.

    scheme = 'explicit' | 'implicit' | 'symmetric'.
    Возвращает словарь с полями 'x', 'v' (решение в последний момент),
    'eps' (относительная погрешность по всей сеточной области),
    'nsteps', 'blown' (признак разрушения счёта), 'frames' -- список
    (t, v) в моменты snapshots, и, если norm_every > 0, историю нормы
    решения 'thist', 'nhist'.
    """
    theta = {'explicit': 0.0, 'implicit': 1.0, 'symmetric': 0.5}[scheme]

    x = np.linspace(pb.A, pb.B, M + 1)
    h = (pb.B - pb.A) / M
    dat = Data(x)
    r = tau * dat.rho / h ** 2

    v = pb.phi(x)
    nsteps = int(np.ceil(T / tau))

    solver = None
    k0 = kM = 0.0
    if theta > 0:
        solver, (k0, kM) = build_matrix(x, h, tau, theta)

    L0, LM = boundary_rows(h)
    lap = np.zeros_like(v)

    num = 0.0   # max |u - v|
    den = 0.0   # max |u|
    if track_error:
        ue = dat.u(0.0)
        num = np.max(np.abs(ue - v))
        den = np.max(np.abs(ue))

    # контрольная оценка: погрешность в фиксированные моменты времени,
    # не зависящие от tau (между слоями решение восстанавливается
    # линейной интерполяцией). Нужна потому, что при больших tau сетка
    # Omega_tau содержит мало слоёв и оценка eps становится оптимистичной.
    cnum = cden = 0.0
    ci = 0
    if ctrl is not None:
        ctrl = np.asarray(ctrl)
        while ci < ctrl.size and ctrl[ci] <= 0.0:
            ue = dat.u(ctrl[ci])
            cnum = max(cnum, np.max(np.abs(ue - v)))
            cden = max(cden, np.max(np.abs(ue)))
            ci += 1

    frames = []
    want = sorted(snapshots) if snapshots else []
    wi = 0
    blown = False
    t_blow = None
    thist = [0.0]
    nhist = [float(np.max(np.abs(v)))]
    ehist = [float(np.max(np.abs(dat.u(0.0) - v)))]

    for n in range(nsteps):
        t = n * tau
        t_new = t + tau
        v_old = v if ctrl is not None else None

        if theta == 0.0:
            lap[1:-1] = r[1:-1] * (v[:-2] - 2 * v[1:-1] + v[2:])
            v = v + lap + tau * dat.f(t)
            v[0] = (dat.ga(t_new) - L0[1] * v[1] - L0[2] * v[2]) / L0[0]
            v[M] = (dat.gb(t_new)
                    - LM[1] * v[M - 1] - LM[2] * v[M - 2]) / LM[0]
        else:
            rhs = v + tau * ((1.0 - theta) * dat.f(t) + theta * dat.f(t_new))
            if theta < 1.0:
                lap[1:-1] = r[1:-1] * (v[:-2] - 2 * v[1:-1] + v[2:])
                rhs += (1.0 - theta) * lap
            rhs[0] = dat.ga(t_new)
            rhs[M] = dat.gb(t_new)
            # та же комбинация строк, что и при сборке матрицы
            rhs[0] -= k0 * rhs[1]
            rhs[M] -= kM * rhs[M - 1]
            v = solver.solve(rhs)

        if norm_every and (n + 1) % norm_every == 0:
            thist.append(t_new)
            nhist.append(float(np.max(np.abs(v))))
            ehist.append(float(np.max(np.abs(dat.u(t_new) - v))))

        if not np.all(np.isfinite(v)) or np.max(np.abs(v)) > blow_at:
            blown = True
            t_blow = t_new
            break

        if track_error:
            ue = dat.u(t_new)
            num = max(num, np.max(np.abs(ue - v)))
            den = max(den, np.max(np.abs(ue)))

        if ctrl is not None:
            while ci < ctrl.size and ctrl[ci] <= t_new + 1e-12:
                s = (ctrl[ci] - t) / tau
                vc = (1.0 - s) * v_old + s * v
                ue = dat.u(ctrl[ci])
                cnum = max(cnum, np.max(np.abs(ue - vc)))
                cden = max(cden, np.max(np.abs(ue)))
                ci += 1

        while wi < len(want) and t_new >= want[wi] - 1e-12:
            frames.append((t_new, v.copy()))
            wi += 1

    out = {'x': x, 'v': v, 'h': h, 'tau': tau, 'nsteps': nsteps,
           'blown': blown, 'frames': frames,
           'thist': np.array(thist), 'nhist': np.array(nhist),
           'ehist': np.array(ehist),
           't_end': t_blow if blown else nsteps * tau,
           'eps': float(num / den) if (track_error and not blown) else np.inf}
    if ctrl is not None:
        out['eps_ctrl'] = float(cnum / cden) if not blown else np.inf
    return out


def tau_stab(h):
    """Условие устойчивости явной схемы: tau <= h^2 / (2 rho_max)."""
    return h ** 2 / (2.0 * pb.RHO_MAX)


def amplification(M, tau):
    """
    Спектральный радиус матрицы перехода явной схемы.

    Оценка tau <= h^2/(2 rho_max) получена для задачи с замороженным
    коэффициентом и без учёта краевых уравнений. Здесь строится матрица
    G перехода со слоя на слой для настоящей задачи:

        строки 1..M-1: G[i,i-1] = r_i, G[i,i] = 1-2 r_i, G[i,i+1] = r_i;
        строка 0 и строка M получаются подстановкой уже пересчитанных
        внутренних значений в краевые уравнения.

    Счёт устойчив тогда и только тогда, когда max|eig(G)| <= 1.
    """
    x = np.linspace(pb.A, pb.B, M + 1)
    h = (pb.B - pb.A) / M
    r = tau * pb.rho(x) / h ** 2

    G = np.zeros((M + 1, M + 1))
    for i in range(1, M):
        G[i, i - 1] = r[i]
        G[i, i] = 1.0 - 2.0 * r[i]
        G[i, i + 1] = r[i]

    L0, LM = boundary_rows(h)
    G[0, :] = -(L0[1] * G[1, :] + L0[2] * G[2, :]) / L0[0]
    G[M, :] = -(LM[1] * G[M - 1, :] + LM[2] * G[M - 2, :]) / LM[0]

    return float(np.max(np.abs(np.linalg.eigvals(G))))
