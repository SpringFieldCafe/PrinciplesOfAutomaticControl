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

exec(compile((SOURCE / 'calculate.py').read_text(encoding='utf-8'), __file__, 'exec'))
out['E7.2']['breakaway_candidates'] = jsroots([4, 18, 26, 10])
out['P7.1']['a']['K_breakaway'] = -float(np.polyval([1, 25, 100, 0], -2.3240812075600177))
out['P7.1']['c']['K_breakaway'] = -float(np.polyval([1, 21, 20, 0], -0.5067535457008205) / (-0.5067535457008205 + 10))
Kcrit = optimize.brentq(lambda x: max(np.roots(add(do, no, x)).real), 1, 2)
out['P7.5']['K1_critical'] = Kcrit
boundx = np.roots([1e-09, -0.0010000008, -2000.0008, 0.4 * ka])
bd = []
for x in boundx:
    if abs(x.imag) < 1e-07 and x.real > 0:
        x = x.real
        k2 = (-1e-06 * x * x + 5.0000008 * x - 0.8) / (2 * ka)
        if k2 > 0:
            bd.append({'K2': k2, 'omega': np.sqrt(x)})
out['P7.11']['K2_imaginary_crossings'] = bd
for name, d, n, label, mark in [('P7_11a', dt, nt, 'Ka, K2 = 10', (ka, 'Ka = 1000')), ('P7_11c', dd, nn, 'K2, Ka = 1000', (10, 'K2 = 10'))]:
    fig, axs = plt.subplots(1, 2, figsize=(10.2, 4.5), layout='constrained')
    kk = logk(-7, 8 if name.endswith('a') else 5)
    locus(axs[0], d, n, kk, f'{name}: complete view', (-2500, 500, -2500, 2500), [mark])
    locus(axs[1], d, n, kk, f'{name}: slow-mode detail', (-1.5, 0.5, -1, 1), [mark])
    for ax in axs:
        ax.set_title(label + ('\nFull view' if ax is axs[0] else '\nNear-origin detail'), fontsize=10)
    fig.savefig(FIGURES / f'{name}.png', dpi=200)
    plt.close(fig)
res = {}
for id, den in [('E7.2', add([1, 6, 13, 10, 0], [1], 6.5)), ('P7.5_inner', qi), ('P7.5_outer', add(do, no, 1)), ('P7.11', add(dt, nt, ka))]:
    rr = np.roots(den)
    res[id] = float(max((abs(np.polyval(den, r)) / np.polyval(abs(den), abs(r)) for r in rr)))
out['normalized_root_residuals'] = res
(BASE / 'numeric_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps({'E7.2_break': out['E7.2']['breakaway_candidates'], 'P7.1': out['P7.1'], 'K1crit': Kcrit, 'K2cross': bd, 'residual': res}, indent=2))
