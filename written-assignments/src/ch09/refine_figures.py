"""Numerical calculations and figures for ch09; paths are script-relative."""
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
import os
H = DATA
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25})
d = np.loadtxt(H / 'P9_5_nyquist.csv', delimiter=',', skiprows=1)
w = d[:, 0]
z = d[:, 1] + 1j * d[:, 2]
fig, axs = plt.subplots(1, 2, figsize=(10, 4.1))
for ax in axs:
    ax.plot(z.real, z.imag, label='Positive frequency')
    ax.plot(z.real, -z.imag, ls='--', label='Negative frequency')
    ax.plot(-1, 0, 'rx', ms=8)
    ax.axhline(0, color='k', lw=0.5)
    ax.axvline(0, color='k', lw=0.5)
    ax.set(xlabel='Re L(j omega)', ylabel='Im L(j omega)')
axs[0].set(title='Full curve', xlim=(min(z.real) - 0.3, max(z.real) + 0.3), ylim=(-max(abs(z.imag)) - 0.4, max(abs(z.imag)) + 0.4))
axs[0].legend(fontsize=8)
axs[1].set(title='Detail near -1', xlim=(-1.15, 0.2), ylim=(-0.6, 0.6))
for ax, target in [(axs[0], 0.4), (axs[1], 1.7)]:
    idx = np.argmin(abs(w - target))
    delta = 35 if ax is axs[0] else 10
    for conjugate in [False, True]:
        zz = np.conj(z) if conjugate else z
        i1, i2 = (idx + delta, idx - delta) if conjugate else (idx - delta, idx + delta)
        ax.annotate('', xy=(zz[i2].real, zz[i2].imag), xytext=(zz[i1].real, zz[i1].imag), arrowprops={'arrowstyle': '->', 'color': 'C1' if conjugate else 'C0', 'lw': 1.8})
fig.tight_layout()
fig.savefig(FIGURES / 'P9_5_nyquist.png', dpi=210)
plt.close(fig)
