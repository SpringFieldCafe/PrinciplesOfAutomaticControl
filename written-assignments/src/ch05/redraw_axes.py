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

from pathlib import Path
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P = DATA
labels = dict(xlabel='$\\mathrm{Re}(s)$ ($\\mathrm{s}^{-1}$)', ylabel='$\\mathrm{Im}(s)$ ($\\mathrm{s}^{-1}$)')
f, ax = plt.subplots(figsize=(5.8, 3.4))
ax.axhline(0, c='k', lw=0.8)
ax.axvline(0, c='k', lw=0.8)
ax.plot([-10], [0], 'x', ms=12, mew=2)
ax.annotate('double pole', (-10, 0), xytext=(-8, 0.6), arrowprops={'arrowstyle': '->'})
ax.set(**labels, xlim=(-15, 2), ylim=(-1.5, 1.5))
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
ax.set(**labels, xlim=(-25, 1), ylim=(-50, 50))
ax.grid(alpha=0.3)
ax.legend()
f.tight_layout()
f.savefig(FIGURES / 'P5.2_region.png', dpi=220)
plt.close(f)
