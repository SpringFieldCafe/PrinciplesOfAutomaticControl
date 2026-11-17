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
import sympy as sp, json
BASE = DATA
a = sp.symbols('a', positive=True, real=True)
s = a * (-1 + 4 * sp.I / 3)
D = s * (s / 1000 + 1) ** 2 * (s * s + 2000 * s + sp.Rational(4, 5))
N = s * s / 50 + 20 * s + sp.Rational(2, 5)
P = sp.factor(sp.im(-D * sp.conjugate(N)))
print('damping intersection equation:', P)
poly = sp.Poly(P / a, a)
inter = []
for v in sp.nroots(poly, maxsteps=300):
    if abs(sp.im(v)) < 1e-10 and sp.re(v) > 0:
        k = complex((-D / N).subs(a, sp.re(v)).evalf())
        if k.real > 0:
            inter.append({'sigma': float(sp.re(v)), 'Ka': k.real, 'imaginary_residual': k.imag})
print(inter)
o = json.loads((BASE / 'numeric_results.json').read_text(encoding='utf-8'))
o['P7.11']['damping_intersection_polynomial'] = str(P)
o['P7.11']['all_positive_damping_intersections'] = inter
o['P7.11']['zeta_boundaries'] = [v['Ka'] for v in inter]
(BASE / 'numeric_results.json').write_text(json.dumps(o, ensure_ascii=False, indent=2), encoding='utf-8')
