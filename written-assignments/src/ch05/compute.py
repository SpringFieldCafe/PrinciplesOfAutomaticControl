"""Numerical calculations and figures for ch05; paths are script-relative."""
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

import json
from pathlib import Path
import numpy as np
import sympy as sp
from scipy import signal, optimize
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
HERE = DATA

def plotlines(name, t, ys, labels, ylabel='$y(t)$', limits=None):
    f, ax = plt.subplots(figsize=(6.8, 3.7))
    for y, l in zip(ys, labels):
        ax.plot(t, y, label=l, lw=1.8)
    ax.set(xlabel='$t$ (s)', ylabel=ylabel, xlim=(t[0], t[-1]))
    ax.grid(alpha=0.3)
    ax.legend()
    f.tight_layout()
    f.savefig(FIGURES / name, dpi=220)
    plt.close(f)
checks = {}
f, ax = plt.subplots(figsize=(5.8, 3.4))
ax.axhline(0, c='k', lw=0.8)
ax.axvline(0, c='k', lw=0.8)
ax.plot([-10], [0], 'x', ms=12, mew=2)
ax.annotate('double pole', (-10, 0), xytext=(-8, 0.6), arrowprops={'arrowstyle': '->'})
ax.set(xlabel='$\\mathrm{Re}(s)$ ($\\mathrm{s}^{-1}$)', ylabel='$\\mathrm{Im}(s)$ ($\\mathrm{s}^{-1}$)', xlim=(-15, 2), ylim=(-1.5, 1.5))
ax.grid(alpha=0.3)
f.tight_layout()
f.savefig(FIGURES / 'E5.6_poles.png', dpi=220)
plt.close(f)
f, ax = plt.subplots(figsize=(6.5, 4))
sig = np.linspace(20 / 3, 25, 300)
lo = np.pi / np.log(10)
hi = np.pi / np.log(5)
for sign in [1, -1]:
    ax.fill_between(-sig, sign * lo * sig, sign * hi * sig, alpha=0.3, color='tab:blue')
ax.axvline(-20 / 3, c='tab:red', ls='--', label='$\\mathrm{Re}(s)=-20/3$')
ax.axhline(0, c='k', lw=0.8)
ax.axvline(0, c='k', lw=0.8)
ax.set(xlabel='$\\mathrm{Re}(s)$ ($\\mathrm{s}^{-1}$)', ylabel='$\\mathrm{Im}(s)$ ($\\mathrm{s}^{-1}$)', xlim=(-25, 1), ylim=(-50, 50))
ax.grid(alpha=0.3)
ax.legend()
f.tight_layout()
f.savefig(FIGURES / 'P5.2_region.png', dpi=220)
plt.close(f)
f, ax = plt.subplots(figsize=(10, 2.5))
ax.set_xlim(0, 13)
ax.set_ylim(-1.5, 1)
ax.axis('off')

def arrow(a, b):
    ax.annotate('', xy=b, xytext=a, arrowprops=dict(arrowstyle='->', lw=1.5))
ax.add_patch(Circle((1.5, 0), 0.17, fill=False))
arrow((0.1, 0), (1.33, 0))
ax.text(0.1, 0.25, '$V_a(s)$')
ax.text(1.25, 0.2, '+')
ax.text(1.15, -0.3, '−')
for xx, ww, label in [(3, 2, '$1/(L_as+R_a)$'), (6, 1, '$K_m$'), (8, 2, '$1/(Js+b)$')]:
    ax.add_patch(Rectangle((xx, -0.35), ww, 0.7, fill=False))
    ax.text(xx + ww / 2, 0, label, ha='center', va='center')
arrow((1.67, 0), (3, 0))
arrow((5, 0), (6, 0))
ax.text(5.1, 0.25, '$I_a$')
arrow((7, 0), (8, 0))
arrow((10, 0), (12.5, 0))
ax.text(11.3, 0.25, '$\\Omega(s)$')
ax.plot([11, 11, 7], [0, -1, -1], 'k')
ax.add_patch(Rectangle((5.5, -1.3), 1.5, 0.6, fill=False))
ax.text(6.25, -1, '$AK_b$', ha='center', va='center')
arrow((5.5, -1), (1.5, -1))
ax.plot([1.5, 1.5], [-1, -0.17], 'k')
f.savefig(FIGURES / 'P5.10_block.png', dpi=220, bbox_inches='tight')
plt.close(f)
t = np.linspace(0, 2, 20001)
ye = 1 - np.exp(-50 * t) / 41 - 40 / 41 * np.exp(-5 * t) * np.cos(5 * t) - 50 / 41 * np.exp(-5 * t) * np.sin(5 * t)
ya = 1 - np.exp(-5 * t) * (np.cos(5 * t) + np.sin(5 * t))
_, yn = signal.step(([2500], [1, 60, 550, 2500]), T=t)
plotlines('E5.15_compare.png', t, [ye, ya], ['Full third order', 'Dominant second order'])
np.savetxt(HERE / 'E5.15_data.csv', np.c_[t, ye, ya], delimiter=',', header='t_s,full,approx', comments='')
checks['E5.15_analytic_scipy_max_error'] = float(np.max(abs(ye - yn)))
checks['E5.15_max_approx_error'] = float(np.max(abs(ye - ya)))
t = np.linspace(0, 20, 20001)
a = 0.51
b = 0.51
wd = np.sqrt(b - a * a / 4)
yst = 1 - np.exp(-a * t / 2) * (np.cos(wd * t) + a / 2 / wd * np.sin(wd * t))
yr = t - 1 + np.exp(-a * t / 2) * (np.cos(wd * t) - 0.745 / wd * np.sin(wd * t))
sys = signal.TransferFunction([b], [1, a, b])
_, ysnum = signal.step(sys, T=t)
_, yrnum, _ = signal.lsim(sys, T=t, U=t)
f, axes = plt.subplots(1, 2, figsize=(10, 3.7))
axes[0].plot(t, yst, label='Step output')
axes[0].axhline(1, ls='--', c='grey')
axes[1].plot(t, t, ls='--', label='Ramp reference')
axes[1].plot(t, yr, label='Ramp output')
axes[1].plot(t, t - yr, label='Ramp error')
for ax in axes:
    ax.set(xlabel='$t$ (s)', ylabel='Response', xlim=(0, 20))
    ax.grid(alpha=0.3)
    ax.legend()
f.tight_layout()
f.savefig(FIGURES / 'AP5.6_responses.png', dpi=220)
plt.close(f)
np.savetxt(HERE / 'AP5.6_data.csv', np.c_[t, yst, yr, t - yr], delimiter=',', header='t_s,step,ramp,ramp_error', comments='')
checks['AP5.6_step_max_error'] = float(np.max(abs(yst - ysnum)))
checks['AP5.6_ramp_max_error'] = float(np.max(abs(yr - yrnum)))

def peak_metrics(K, detail=False):
    den = [1, 6, 15 + K, 50 + 2 * K]
    r, ps, _ = signal.residue([K, 2 * K], np.polymul(den, [1, 0]))
    rr, pp, _ = signal.residue([K, 2 * K], den)
    final = 2 * K / (50 + 2 * K)
    horizon = 30 / min(-pp.real)
    freq = max(abs(pp.imag))
    dt = min(0.02, np.pi / max(freq, 1) / 30)
    tt = np.arange(0, horizon + dt, dt)
    v = (rr[:, None] * np.exp(pp[:, None] * tt)).sum(axis=0).real

    def dy(t):
        return np.real(np.sum(rr * np.exp(pp * t)))
    roots = [optimize.brentq(dy, tt[i], tt[i + 1]) for i in np.where(v[:-1] * v[1:] < 0)[0]]

    def yy(t):
        return np.real(np.sum(r * np.exp(ps * t)))
    peak = max([(yy(t), t) for t in roots] + [(final, horizon)])
    po = 100 * (peak[0] / final - 1)
    return (po, peak, final) if detail else po
ks = np.logspace(-5, 5, 301)
pos = np.array([peak_metrics(k) for k in ks])
cand = np.where((pos[1:-1] < pos[:-2]) & (pos[1:-1] < pos[2:]))[0] + 1
opts = [optimize.minimize_scalar(peak_metrics, bounds=(ks[i - 1], ks[i + 1]), method='bounded', options={'xatol': 1e-08}) for i in cand]
best = min(opts, key=lambda r: r.fun)
Kopt = best.x
po, pk, final = peak_metrics(Kopt, True)
t = np.linspace(0, 3, 30001)
_, y = signal.step(([Kopt, 2 * Kopt], [1, 6, 15 + Kopt, 50 + 2 * Kopt]), T=t)
f, axes = plt.subplots(1, 2, figsize=(10, 3.7))
axes[0].semilogx(ks, pos)
axes[0].plot([Kopt], [po], 'ro')
axes[0].set(xlabel='$K$', ylabel='Overshoot (%)')
axes[1].plot(t, y, label='Optimal step response')
axes[1].axhline(final, c='grey', ls='--', label='Steady state')
axes[1].plot([pk[1]], [pk[0]], 'ro')
axes[1].set(xlabel='$t$ (s)', ylabel='$y(t)$')
axes[1].legend()
for ax in axes:
    ax.grid(alpha=0.3)
f.tight_layout()
f.savefig(FIGURES / 'P5.15_optimum.png', dpi=220)
plt.close(f)
np.savetxt(HERE / 'P5.15_scan.csv', np.c_[ks, pos], delimiter=',', header='K,PO_percent', comments='')
np.savetxt(HERE / 'P5.15_step.csv', np.c_[t, y], delimiter=',', header='t_s,y', comments='')
checks['P5.15'] = {'K': Kopt, 'PO_percent': po, 'peak': pk, 'steady': final, 'scipy_sample_peak': float(y.max()), 'candidate_count': len(opts), 'K_samples': {str(k): peak_metrics(k) for k in [20, 50, 100, 200]}}

def perf(q):
    w2 = 1 + q * q
    num = [5 * w2]
    den = np.polymul([1, 5], [1, 2, w2])
    r, pp, _ = signal.residue(num, np.polymul(den, [1, 0]))
    rr, ppd, _ = signal.residue(num, den)

    def yy(t):
        return np.real(np.sum(r * np.exp(pp * t)))

    def dy(t):
        return np.real(np.sum(rr * np.exp(ppd * t)))
    tp = optimize.brentq(dy, 0.1, 4 / q)
    tt = np.linspace(0.001, 15, 10001)
    v = (r[:, None] * np.exp(pp[:, None] * tt)).sum(axis=0).real
    ii = np.where(abs(v - 1) > 0.02)[0][-1]
    target = 1 + (0.02 if v[ii] > 1 else -0.02)
    ts = optimize.brentq(lambda u: yy(u) - target, tt[ii], tt[ii + 1])
    return (yy(tp) - 1, tp, ts)
q = optimize.brentq(lambda q: perf(q)[0] - 0.2, 1.9, 2.4)
mp, tp, us = perf(q)
sig = us / 0.6
w2 = sig * sig * (1 + q * q)
r3 = 5 * sig
checks['P5.2_exact'] = {'q': q, 'sigma': sig, 'omega_d': q * sig, 'r3': r3, 'wn_squared': w2, 'G_numerator': r3 * w2, 'G_den_s2_coefficient': r3 + 2 * sig, 'G_den_s_coefficient': w2 + 2 * sig * r3, 'PO_percent': mp * 100, 'Ts': us / sig}
x = sp.symbols('s')
GH = (x ** 3 + 7 * x * x + 24 * x + 24) / ((x + 1) * (x + 2) * (x + 3) * (x + 4))
GL = (sp.Rational(7, 24) * x + 1) / (sp.Rational(115, 288) * x * x + sp.Rational(11, 8) * x + 1)
checks['P5.13_residual_series'] = str(sp.series(GH - GL, x, 0, 4))
checks['P5.10_A_bound'] = np.pi ** 2 / np.log(0.15) ** 2
(HERE / 'numeric_results.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding='utf-8')
print(checks)
