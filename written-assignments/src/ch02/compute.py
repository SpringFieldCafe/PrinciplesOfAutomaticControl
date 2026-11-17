"""Numerical calculations and figures for ch02; paths are script-relative."""
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

import json, math
from pathlib import Path
import sympy as s
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
HERE = DATA
fig, ax = plt.subplots(figsize=(10, 2.7))
ax.set_xlim(0, 12)
ax.set_ylim(-1.7, 1.1)
ax.axis('off')

def arrow(a, b):
    ax.annotate('', xy=b, xytext=a, arrowprops=dict(arrowstyle='->', lw=1.5))
ax.add_patch(Circle((1.6, 0), 0.18, fill=False))
ax.text(1.4, 0.25, '+', fontsize=12)
ax.text(1.33, -0.34, '−', fontsize=12)
arrow((0.1, 0), (1.4, 0))
ax.text(0.1, 0.22, '$R(s)$')
for x, w, label in [(3, 1.3, '$8$'), (5.4, 1.6, '$0.5/s$'), (8, 2.4, '$1/(s^2+2s+4)$')]:
    ax.add_patch(Rectangle((x, -0.4), w, 0.8, fill=False, lw=1.4))
    ax.text(x + w / 2, 0, label, ha='center', va='center', fontsize=13)
arrow((1.8, 0), (3, 0))
ax.text(2, 0.22, '$V_1(s)$')
arrow((4.3, 0), (5.4, 0))
ax.text(4.5, 0.22, '$V_2(s)$')
arrow((7, 0), (8, 0))
ax.text(7.2, 0.22, '$\\Theta(s)$')
arrow((10.4, 0), (11.8, 0))
ax.text(11, 0.22, '$P(s)$')
ax.plot([11, 11, 1.6, 1.6], [0, -1.15, -1.15, -0.18], 'k-', lw=1.4)
ax.annotate('', xy=(5, -1.15), xytext=(7, -1.15), arrowprops=dict(arrowstyle='->', lw=1.4))
fig.savefig(FIGURES / 'E2.16_block.png', dpi=220, bbox_inches='tight')
plt.close(fig)
R = 10000 * math.exp(-2)
ss = s.symbols('s')
M, m, b, k = s.symbols('M m b k')
A, B, C, D, H1, H2 = s.symbols('A B C D H1 H2')
matrix = s.Matrix([[1 + A * H2, B * H1], [C * H2, 1 + D * H1]])
delta = matrix.det()
out = matrix.inv() * s.Matrix([A - B, C - D])
checks = {'P2.32_symbolic_residual': str(s.simplify(out[1] - (C - D - (A * D - B * C) * H2) / delta)), 'P2.26_determinant_residual': str(s.expand((M * ss ** 2 + b * ss + k) * (m * ss ** 2 + b * ss + k) - (b * ss + k) ** 2 - ss ** 2 * (M * m * ss ** 2 + (M + m) * b * ss + (M + m) * k))), 'E2.2_Rstar': R, 'E2.2_derivative': -0.1 * R, 'E2.2_finite_difference': (10000 * math.exp(-0.1 * (20 + 0.0001)) - 10000 * math.exp(-0.1 * (20 - 0.0001))) / 0.0002}
(HERE / 'numeric_results.json').write_text(json.dumps(checks, ensure_ascii=False, indent=2), encoding='utf-8')
print(checks)
