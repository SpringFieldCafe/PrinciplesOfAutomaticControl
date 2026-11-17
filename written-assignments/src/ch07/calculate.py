"""Numerical calculations and figures for ch07; paths are script-relative."""
from pathlib import Path
import os
SOURCE = Path(__file__).resolve().parent
ROOT = SOURCE.parents[1]
DATA = ROOT / 'data' / SOURCE.name
FIGURES = ROOT / 'figures' / SOURCE.name
for directory in (DATA, FIGURES):
    directory.mkdir(parents=True, exist_ok=True)
os.environ['MPLCONFIGDIR'] = str(ROOT / '.local' / 'cache' / 'matplotlib')
Path(os.environ['MPLCONFIGDIR']).mkdir(parents=True, exist_ok=True)

from pathlib import Path
import os, json, itertools
BASE = DATA
import numpy as np
import sympy as sp
from scipy import optimize, signal
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
s = sp.symbols('s')

def jsroots(x):
    return [[float(z.real), float(z.imag)] for z in np.roots(x)]

def add(d, n, k):
    return np.polyadd(d, k * np.asarray(n))

def locus(ax, d, n, ks, title, lim, marks=None):
    nroot = len(d) - 1
    rr = []
    prev = None
    for k in ks:
        v = np.roots(add(d, n, k))
        if prev is not None:
            ix = linear_sum_assignment(abs(prev[:, None] - v[None, :]))[1]
            v = v[ix]
        rr.append(v)
        prev = v
    rr = np.asarray(rr)
    for j in range(nroot):
        line, = ax.plot(rr[:, j].real, rr[:, j].imag, lw=1.2)
        v = rr[:, j]
        inside = np.where((v.real > lim[0]) & (v.real < lim[1]) & (v.imag > lim[2]) & (v.imag < lim[3]))[0]
        if len(inside) > 30:
            ii = inside[int(len(inside) * 0.7)]
            jj = min(ii + 18, len(v) - 1)
            ax.annotate('', xy=(v[jj].real, v[jj].imag), xytext=(v[ii].real, v[ii].imag), arrowprops={'arrowstyle': '->', 'color': line.get_color(), 'lw': 1.3})
    poles = np.roots(d)
    zeros = np.roots(np.trim_zeros(n, 'f'))
    ax.scatter(poles.real, poles.imag, marker='x', color='black', s=40, label='Open-loop poles', zorder=5)
    if len(zeros):
        ax.scatter(zeros.real, zeros.imag, marker='o', facecolors='none', edgecolors='black', s=40, label='Zeros', zorder=5)
    nd = len(poles) - len(zeros)
    if nd > 0:
        center = (sum(poles) - sum(zeros)).real / nd
        length = 3 * max(lim[1] - lim[0], lim[3] - lim[2])
        tt = np.linspace(0, length, 50)
        for j in range(nd):
            th = np.pi * (2 * j + 1) / nd
            ax.plot(center + tt * np.cos(th), tt * np.sin(th), color='gray', ls='--', lw=0.7, alpha=0.65, label='Asymptotes' if j == 0 else None)
    if marks is not None:
        for k, lab in marks:
            v = np.roots(add(d, n, k))
            ax.scatter(v.real, v.imag, s=26, label=lab, zorder=6)
    ax.axhline(0, color='gray', lw=0.6)
    ax.axvline(0, color='gray', lw=0.6)
    ax.grid(True, alpha=0.25)
    ax.set(xlabel='Re(s) (1/s)', ylabel='Im(s) (1/s)', title=title, xlim=lim[:2], ylim=lim[2:])
    ax.legend(fontsize=7, loc='best')
    np.savetxt(BASE / (title.split(':')[0].replace('.', '_').replace(' ', '_') + '_locus.csv'), np.c_[ks, rr.real, rr.imag], delimiter=',', header='K,' + ','.join((f're{i}' for i in range(nroot))) + ',' + ','.join((f'im{i}' for i in range(nroot))), comments='')

def one(name, d, n, ks, title, lim, marks=None):
    fig, ax = plt.subplots(figsize=(7, 4.6), layout='constrained')
    locus(ax, d, n, ks, title, lim, marks)
    fig.savefig(FIGURES / f'{name}.png', dpi=200)
    plt.close(fig)
logk = lambda lo, hi: np.r_[0, np.logspace(lo, hi, 2500)]
out = {}
K = (np.sqrt(5) - 1) / 4
p = -(1 + 2 * K) / (1 + K)
out['E7.1'] = {'K': K, 'repeated_pole': p, 'ts_2percent_pole_estimate': 4 / -p, 'roots': jsroots(add([1, 2, 2], [1, 4, 0], K))}
one('E7_1', [1, 2, 2], [1, 4, 0], logk(-5, 5), 'E7.1: root locus, K >= 0', (-4.5, 0.5, -1.6, 1.6), [(K, 'Repeated root')])
d = [1, 6, 13, 10, 0]
r = np.roots(add(d, [1], 6.5))
pair = sorted([z for z in r if z.imag > 0], key=lambda z: z.real, reverse=True)[0]
sig = -pair.real
wd = pair.imag
out['E7.2'] = {'roots': jsroots(add(d, [1], 6.5)), 'zeta': sig / abs(pair), 'ts': 4 / sig, 'overshoot_percent': 100 * np.exp(-np.pi * sig / wd)}
one('E7_2', d, [1], logk(-4, 4), 'E7.2: root locus', (-5, 2, -4, 4), [(6.5, 'K = 6.5')])
d = [1, 50, 500, 1000]
st = np.roots(np.polyder(d))
out['E7.5'] = {'open_poles': jsroots(d), 'stationary': [[float(x), float(-np.polyval(d, x))] for x in st], 'centroid': -50 / 3}
out['E7.9'] = {'centroid': -2 / 3, 'upper_departure_deg': -np.degrees(np.arctan(0.5)), 'Kcrit': 10, 'roots_critical': jsroots([1, 2, 5, 10])}
one('E7_9', [1, 2, 5, 0], [1], logk(-5, 4), 'E7.9: root locus', (-8, 5, -8, 8), [(10, 'K = 10')])
cases = [([1, 25, 100, 0], [1], (-30, 10, -25, 25)), ([1, 4, 6, 4], [1], (-7, 4, -6, 6)), ([1, 21, 20, 0], [1, 10], (-25, 3, -25, 25)), ([1, 1, 0, 0], [1, 4, 8], (-8, 3, -5, 5))]
fig, axs = plt.subplots(2, 2, figsize=(10, 7.6), layout='constrained')
for i, (d, n, lim) in enumerate(cases):
    locus(axs.flat[i], d, n, logk(-5, 6), f'P7.1({chr(97 + i)}): root locus', lim)
fig.savefig(FIGURES / 'P7_1.png', dpi=200)
plt.close(fig)
out['P7.1'] = {}
for i, (d, n, lim) in enumerate(cases):
    poly = sp.Poly(sum((float(x) * s ** (len(d) - 1 - j) for j, x in enumerate(d))), s).as_expr()
    num = sum((float(x) * s ** (len(n) - 1 - j) for j, x in enumerate(n)))
    der = sp.together(sp.diff(-poly / num, s))
    c = sp.Poly(sp.fraction(der)[0], s)
    out['P7.1'][chr(97 + i)] = {'stationary': [[float(sp.re(x)), float(sp.im(x))] for x in sp.nroots(c)], 'poles': jsroots(d), 'zeros': jsroots(n)}
a4 = 0.155
a3 = 1 + 0.155 * 0.707
a2 = 0.707 + 0.155 * 0.25
a1 = 0.25
kc = (a3 * a2 * a1 - a4 * a1 * a1) / (a3 * a3 * 0.25)
d = [a4, a3, a2, a1, 0]
out['P7.4'] = {'denominator': d, 'ka_critical': kc, 'omega_critical': np.sqrt(a1 / a3), 'roots': jsroots(add(d, [0.25], kc))}
one('P7_4', d, [0.25], logk(-5, 5), 'P7.4: antenna root locus', (-8, 3, -5, 5), [(kc, f'ka = {kc:.5f}')])
dg = np.polymul([1, 0.4], [1, -0.36, 0.16])
ng = [25, 0.75]
di = np.polymul(dg, [1, 9])
ni = np.polymul(ng, [1, 1])
z = 0.707
ratio = np.sqrt(1 - z * z) / z

def kr(a):
    v = -a + 1j * a * ratio
    return -np.polyval(di, v) / np.polyval(ni, v)
grid = np.logspace(-5, 4, 30000)
inter = []
for a, b in zip(grid[:-1], grid[1:]):
    if kr(a).imag * kr(b).imag < 0:
        rt = optimize.brentq(lambda x: kr(x).imag, a, b)
        kv = kr(rt).real
        if kv > 0:
            inter.append((rt, kv, np.roots(add(di, ni, kv))))
out['P7.5'] = {'zeta_intersections': [{'a': a, 'K2': k, 'roots': [[v.real, v.imag] for v in rr]} for a, k, rr in inter]}
sel = [x for x in inter if max(x[2].real) < 0][0]
K2 = sel[1]
qi = add(di, ni, K2)
out['P7.5'].update({'K2': K2, 'inner_coefficients': qi.tolist(), 'disturbance_y_inner': 6.75 / qi[-1]})
one('P7_5a', di, ni, logk(-6, 5), 'P7.5(a): automatic stabilization', (-14, 3, -12, 12), [(K2, f'K2 = {K2:.5f}')])
do = np.polymul([1, 12, 1], qi)
no = np.polymul(ng, [1, 9])
test = []
for k in [0.01, 0.1, 0.5, 1, 2, 5, 10]:
    rr = np.roots(add(do, no, k))
    test.append({'K1': k, 'max_real': float(max(rr.real))})
K1 = 1
out['P7.5'].update({'K1_tests': test, 'K1_chosen': K1, 'outer_denominator': do.tolist(), 'outer_numerator': no.tolist(), 'outer_roots': jsroots(add(do, no, K1)), 'disturbance_y_outer': 6.75 / (qi[-1] + 6.75 * K1)})
one('P7_5c', do, no, logk(-6, 6), 'P7.5(c): pilot-loop root locus', (-15, 4, -10, 10), [(K1, 'K1 = 1')])
dt = np.polymul(np.polymul([0.001, 1], [0.001, 1]), [1, 2000, 0.8])
dt = np.r_[dt, 0]
nt = [0.02, 20, 0.4]

def minz(k):
    r = np.roots(add(dt, nt, k))
    return float(min(-r.real / abs(r)))
gg = np.logspace(-7, 8, 18000)
intervals = []
vv = [minz(k) - 0.6 for k in gg]
for a, b, fa, fb in zip(gg[:-1], gg[1:], vv[:-1], vv[1:]):
    if fa * fb < 0:
        intervals.append(optimize.brentq(lambda x: minz(x) - 0.6, a, b, xtol=1e-09))
out['P7.11'] = {'D0': dt.tolist(), 'NKa': nt, 'zeta_boundaries': intervals, 'sample_damping': [[float(k), minz(k)] for k in [0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]]}
one('P7_11a', dt, nt, logk(-6, 8), 'P7.11(a): Ka root locus, K2 = 10', (-2500, 500, -2500, 2500))
ka = 1000
dd = np.polyadd(dt, [0.4 * ka])
nn = [0.002 * ka, 2 * ka, 0]
out['P7.11'].update({'Ka_chosen': ka, 'roots_chosen': jsroots(add(dt, nt, ka)), 'min_zeta_chosen': minz(ka), 'fixed_Ka_D': dd.tolist(), 'NK2': nn})
one('P7_11c', dd, nn, logk(-6, 5), 'P7.11(c): K2 root locus, Ka = 1000', (-2500, 500, -2500, 2500), [(10, 'K2 = 10')])
(BASE / 'numeric_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(out, indent=2))
