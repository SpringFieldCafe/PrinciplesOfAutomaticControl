"""Numerical calculations and figures for ch10; paths are script-relative."""
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
d = np.loadtxt(H / 'P10_43_overshoot.csv', delimiter=',', skiprows=1)
kl = 35 - 4 * np.sqrt(76)
ku = 35 + 4 * np.sqrt(76)
low = d[:, 0] < kl
fig, ax = plt.subplots(1, 2, figsize=(9, 3.9))
ax[0].plot(d[low, 0], d[low, 1])
ax[1].plot(d[~low, 0], d[~low, 1])
ax[0].set(xlim=(0, kl), xlabel='Gain K', ylabel='Overshoot (%)', title='Low-gain stable interval')
ax[1].set(xlim=(ku, 100), xlabel='Gain K', ylabel='Overshoot (%)', title='High-gain stable interval')
fig.suptitle('The middle interval 0.128808 < K < 69.871192 is unstable', fontsize=11)
fig.tight_layout()
fig.savefig(FIGURES / 'P10_43_overshoot.png', dpi=210)
plt.close(fig)
d = np.loadtxt(H / 'P10_19_step.csv', delimiter=',', skiprows=1)
fig, ax = plt.subplots(figsize=(7.1, 4))
ax.plot(d[:, 0], d[:, 1], label='Uncompensated, K=1')
ax.plot(d[:, 0], d[:, 2], label='Lag compensated, K=1, b=1/15')
ax.axhline(1, color='k', lw=0.8, label='Reference')
ax.axhline(10 / 11, color='C1', ls=':', lw=0.8)
ax.axhline(0.25, color='C0', ls=':', lw=0.8)
ax.set(xlim=(0, 18), xlabel='Time (s)', ylabel='Output')
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIGURES / 'P10_19_step.png', dpi=210)
plt.close(fig)
d = np.loadtxt(H / 'P10_1_rootlocus.csv', delimiter=',', skiprows=1)
fig, ax = plt.subplots(figsize=(6.8, 4))
for j in range(3):
    ax.scatter(d[:, 1 + j], d[:, 4 + j], s=1, color='C0')
ax.plot([-2.4, -2.4, -12], [3.2, -3.2, 0], 'rx', ms=8, label='Selected closed-loop poles')
ax.set(xlim=(-20, 2), ylim=(-12, 12), xlabel='Real part (1/s)', ylabel='Imaginary part (rad/s)')
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES / 'P10_1_rootlocus.png', dpi=210)
plt.close(fig)
