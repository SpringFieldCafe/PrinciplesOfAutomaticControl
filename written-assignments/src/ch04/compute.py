"""Numerical calculations and figures for ch04; paths are script-relative."""
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
import sympy as s
from scipy import signal
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = DATA
t = np.linspace(0, 10, 2001)
rt = np.sqrt(1905)
alpha = (55 - rt) / 14
beta = (55 + rt) / 14
er = -7 / 8 + (15 / 16 + 393 / (16 * rt)) * np.exp(-alpha * t) + (15 / 16 - 393 / (16 * rt)) * np.exp(-beta * t)
ed = 7 / 4 * (1 - (beta * np.exp(-alpha * t) - alpha * np.exp(-beta * t)) / (beta - alpha))
_, num_er = signal.step(signal.TransferFunction([7, 28, -35], [7, 55, 40, 0]), T=t)
_, num_ed = signal.step(signal.TransferFunction([70], [7, 55, 40]), T=t)
_, num_er = signal.step(signal.TransferFunction([7, 28, -35], [7, 55, 40]), T=t)
for arr, name, limit, title in [(er, 'AP4.3_tracking.png', -0.875, 'Tracking error'), (ed, 'AP4.3_disturbance.png', 1.75, 'Disturbance error')]:
    fig0, ax = plt.subplots(figsize=(6.8, 3.6))
    ax.plot(t, arr, lw=2, label='$e(t)$')
    ax.axhline(limit, ls='--', c='grey', label='Steady state')
    ax.set(xlabel='$t$ (s)', ylabel='$e(t)$', xlim=(0, 10), title=title)
    ax.grid(alpha=0.3)
    ax.legend()
    fig0.tight_layout()
    fig0.savefig(FIGURES / name, dpi=220)
    plt.close(fig0)
np.savetxt(HERE / 'AP4.3_data.csv', np.c_[t, er, ed], delimiter=',', header='t_s,tracking_error,disturbance_error', comments='')
x, K, b = s.symbols('s K b')
checks = {'AP4.3_tracking_max_error': float(np.max(abs(er - num_er))), 'AP4.3_disturbance_max_error': float(np.max(abs(ed - num_ed))), 'AP4.3_poles': [-float(alpha), -float(beta)], 'P4.8_stable_poles_K1': np.roots([0.1, 3.5, 45, 400]).tolist(), 'P4.8_unstable_poles_K7': [(float(z.real), float(z.imag)) for z in np.roots([0.1, 3.5, 45, 1600])], 'P4.17_wd': float(np.sqrt(2993.75))}
checks['P4.8_stable_poles_K1'] = [(float(z.real), float(z.imag)) for z in np.roots([0.1, 3.5, 45, 400])]
(HERE / 'numeric_results.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding='utf-8')
print(checks)
