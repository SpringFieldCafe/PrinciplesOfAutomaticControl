"""Numerical calculations and figures for ch08; paths are script-relative."""
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
import json, numpy as np
from scipy.optimize import brentq
BASE = DATA
j = json.loads((BASE / 'numeric_results.json').read_text(encoding='utf-8'))
z = j['P8.24']['zeta']
wn = 0.9
wd = wn * np.sqrt(1 - z * z)
error = lambda t: np.exp(-z * wn * t) / np.sqrt(1 - z * z) * np.sin(wd * t + np.arccos(z))
t = np.linspace(0, 100, 100000)
v = abs(error(t)) - 0.02
ind = np.where(v > 0)[0][-1]
ts = brentq(lambda t: abs(error(t)) - 0.02, t[ind], t[ind + 1])
j['P8.24']['exact_2percent_model_ts_s'] = ts
j['E8.5'] = {'K': 8, 'a_seconds': 0.25, 'b_seconds': 1 / 24}
(BASE / 'numeric_results.json').write_text(json.dumps(j, ensure_ascii=False, indent=2), encoding='utf-8')
print(ts)
