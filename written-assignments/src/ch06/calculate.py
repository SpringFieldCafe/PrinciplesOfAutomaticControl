"""Numerical calculations and figures for ch06; paths are script-relative."""
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
import os, json
BASE = DATA
import numpy as np
import sympy as sp
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
s = sp.symbols('s')

def roots(c):
    return [[float(z.real), float(z.imag)] for z in np.roots(c)]
out = {}
z = 0.707
kv = 18000 / np.pi
out['E6.6'] = {'K_zeta': float(np.roots([1, -(2 + 8 * z * z), 1])[0]), 'imaginary_K': 1}
out['E6.11'] = {'roots': roots([1, 6, 2, 1, 3]), 'formal_DC_error': -7 / 3}
q = s ** 6 + 9 * s ** 5 + sp.Rational(125, 4) * s ** 4 + sp.Rational(245, 4) * s ** 3 + sp.Rational(271, 4) * s ** 2 + sp.Rational(59, 4) * s + 15
out['E6.15'] = {'factorization': str(sp.factor(q)), 'roots': roots([1, 9, 31.25, 61.25, 67.75, 14.75, 15])}
out['P6.4'] = {'K_critical': 600, 'omega_critical': float(np.sqrt(800)), 'K_half_roots': roots([1, 30, 500, 12000]), 'shift_1': [1, 27, 443, 11529], 'shift_2': [1, 24, 392, 11112]}
out['P6.7'] = {'Kv': kv, 'roots': roots([1, 101, 100 + 10 * kv, 100 * kv]), 'steady_error_rad': 100 / kv}
out['AP6.3'] = {'alpha_min': float((np.sqrt(5) - 1) / 2), 'selected_roots': roots([1, 2, 1, 1])}
gains = {'K1': 1, 'K2': 1, 'K3': 20, 'K4': 300, 'K5': 46}
out['AP6.5'] = {'a_gains': gains, 'a_denominator': [2, 64, 42, 180], 'a_roots': roots([2, 64, 42, 180]), 'a_DC': 30 / 180, 'b_feasible_positive_gains': False}
t = np.linspace(0, 100, 16001)
tt, y = signal.step(signal.TransferFunction([30], [2, 64, 42, 180]), T=t)
np.savetxt(BASE / 'AP6_5_step.csv', np.c_[tt, y], delimiter=',', header='time_s,y', comments='')
fig, ax = plt.subplots(figsize=(7, 3.8), layout='constrained')
ax.plot(tt, y, label='K1=K2=1, K3=20, K4=300, K5=46')
ax.axhline(1 / 6, color='gray', ls='--', label='Final value = 1/6')
ax.set(xlabel='Time (s)', ylabel='Output y(t)', title='AP6.5(a): stable unit-step response')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)
fig.savefig(FIGURES / 'AP6_5_step.png', dpi=200)
plt.close(fig)
kd = np.linspace(0, 100, 301)
kp = 210 + 9 * kd
np.savetxt(BASE / 'AP6_6_boundary.csv', np.c_[kd, kp], delimiter=',', header='Kd,Kp_boundary', comments='')
fig, ax = plt.subplots(figsize=(6.7, 4), layout='constrained')
ax.plot(kd, kp, label='Kp = 210 + 9 Kd')
ax.fill_between(kd, 0, kp, alpha=0.2, label='Stable: 0 < Kp < 210 + 9 Kd')
ax.set(xlabel='Kd', ylabel='Kp', xlim=(0, 100), ylim=(0, 1200), title='AP6.6: positive-gain stability region')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)
fig.savefig(FIGURES / 'AP6_6_region.png', dpi=200)
plt.close(fig)
(BASE / 'numeric_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(out, indent=2))
