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
import os, json
HERE = DATA
import numpy as np
from scipy import signal, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'savefig.dpi': 210})
RESULT = {}

def tf(n, d):
    return lambda w: np.polyval(n, 1j * np.asarray(w)) / np.polyval(d, 1j * np.asarray(w))

def rootscan(fun, lo=1e-05, hi=100000.0, n=40000):
    w = np.geomspace(lo, hi, n)
    y = fun(w)
    inds = np.where(y[:-1] * y[1:] < 0)[0]
    return [float(optimize.brentq(fun, w[i], w[i + 1], xtol=1e-12)) for i in inds]

def margins(fun, lo=1e-05, hi=100000.0):
    w = np.geomspace(lo, hi, 40000)
    ph = np.unwrap(np.angle(fun(w))) * 180 / np.pi
    wc = rootscan(lambda x: abs(fun(x)) - 1, lo, hi)
    pc = [x for x in rootscan(lambda x: np.imag(fun(x)), lo, hi) if np.real(fun(x)) < 0]

    def exact_phase(x):
        direct = float(np.angle(fun(x), deg=True))
        approx = float(np.interp(np.log(x), np.log(w), ph))
        return direct + 360 * round((approx - direct) / 360)
    return {'gain_crossovers': [{'w': x, 'PM_deg': float(180 + exact_phase(x))} for x in wc], 'phase_crossovers': [{'w': x, 'GM': float(1 / abs(fun(x))), 'GM_dB': float(-20 * np.log10(abs(fun(x))))} for x in pc]}

def closed_metrics(fun):
    w = np.geomspace(1e-05, 100000.0, 50000)
    a = abs(fun(w) / (1 + fun(w)))
    i = np.argmax(a)
    r = optimize.minimize_scalar(lambda x: -abs(fun(np.exp(x)) / (1 + fun(np.exp(x)))), bounds=(np.log(w[max(0, i - 2)]), np.log(w[min(len(w) - 1, i + 2)])), method='bounded')
    peak = -r.fun
    wr = np.exp(r.x)
    bw = rootscan(lambda x: abs(fun(x) / (1 + fun(x))) - 1 / np.sqrt(2))[0]
    return {'Mr': float(peak), 'Mr_dB': float(20 * np.log10(peak)), 'wr': float(wr), 'bandwidth': float(bw)}

def bode(name, funs, labels, w, marks=None):
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.7), sharex=True)
    cols = []
    for f, label in zip(funs, labels):
        z = f(w)
        mag = 20 * np.log10(abs(z))
        phase = np.unwrap(np.angle(z)) * 180 / np.pi
        ax[0].semilogx(w, mag, label=label)
        ax[1].semilogx(w, phase, label=label)
        cols += [mag, phase]
    ax[0].axhline(0, color='k', lw=0.7)
    ax[1].axhline(-180, color='k', lw=0.7)
    ax[0].set_ylabel('Magnitude (dB)')
    ax[1].set_ylabel('Phase (deg)')
    ax[1].set_xlabel('Frequency (rad/s)')
    ax[0].legend(loc='best')
    if marks:
        for wc in marks:
            for a in ax:
                a.axvline(wc, color='grey', ls=':', lw=0.9)
    fig.tight_layout()
    fig.savefig(FIGURES / (name + '.png'))
    plt.close(fig)
    np.savetxt(HERE / (name + '.csv'), np.column_stack([w, *cols]), delimiter=',', header='omega,' + ','.join((f'{label}_{v}' for label in labels for v in ['dB', 'deg'])), comments='')

def step(name, n, d, t):
    t, y = signal.step((n, d), T=t)
    err = abs(y - 1)
    bad = np.where(err > 0.02)[0]
    ts = float(t[bad[-1] + 1]) if len(bad) and bad[-1] + 1 < len(t) else None
    r = {'overshoot_pct': float((max(y) - 1) * 100), 'peak_time': float(t[np.argmax(y)]), 'settling_2pct': ts, 'poles': [[float(x.real), float(x.imag)] for x in np.roots(d)]}
    fig, ax = plt.subplots(figsize=(7.1, 3.8))
    ax.plot(t, y, label='Step response')
    ax.axhline(1, color='k', lw=0.8)
    ax.axhspan(0.98, 1.02, color='grey', alpha=0.14, label='2% band')
    ax.set(xlabel='Time (s)', ylabel='Output')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / (name + '.png'))
    plt.close(fig)
    np.savetxt(HERE / (name + '.csv'), np.c_[t, y], delimiter=',', header='time_s,output', comments='')
    return r

def run():
    f = tf([0.2, 2], np.polymul([1, 0], np.polymul([5, 1], [1 / 81, 1 / 9, 1])))
    RESULT['E9.1'] = margins(f)
    bode('E9_1_bode', [f], ['L'], np.geomspace(0.01, 1000, 4000))
    wn = np.sqrt(150)
    ze = 5 / (2 * wn)
    Mr = 1 / (2 * ze * np.sqrt(1 - ze * ze))
    wr = wn * np.sqrt(1 - 2 * ze * ze)
    bw = wn * np.sqrt(1 - 2 * ze * ze + np.sqrt(2 - 4 * ze * ze + 4 * ze ** 4))
    Mp = 100 * np.exp(-np.pi * ze / np.sqrt(1 - ze * ze))
    RESULT['E9.13'] = {'wn': wn, 'zeta': ze, 'Mr': Mr, 'Mr_dB': 20 * np.log10(Mr), 'wr': wr, 'bandwidth': bw, 'overshoot_pct': Mp, **closed_metrics(tf([150], [1, 5, 0]))}
    base = tf([1, 2], [1, 2, 15, 0])
    wc = rootscan(lambda w: np.unwrap(np.angle(base(np.atleast_1d(w))))[0] + 5 * np.pi / 6 if np.isscalar(w) else np.angle(base(w)) + 5 * np.pi / 6, lo=3.9, hi=20)[0]
    K = 1 / abs(base(wc))
    f = lambda w: K * base(w)
    RESULT['E9.17'] = {'K': K, 'wc': wc, 'Kv': 2 * K / 15, 'ramp_error': 15 / (2 * K), 'margins': margins(f), 'poles': [[x.real, x.imag] for x in np.roots([1, 2, 15 + K, 2 * K])]}
    bode('E9_17_bode', [base, f], ['K=1', f'K={K:.5g}'], np.geomspace(0.01, 1000, 4000), [wc])
    RESULT['P9.3'] = {'roots': [[x.real, x.imag] for x in np.roots([1, 11, 56, 96])], 'shifted_coefficients': [1, 8, 37, 50], 'shifted_routh_first_column': [1, 8, 30.75, 50]}
    den = np.polymul(np.polymul([1, 1], [3, 1]), [0.4, 1])
    f = tf([10], den)
    RESULT['P9.5'] = {'K': 10, 'K_steady_lower': 9, 'K_stability_upper': 4.6 * 4.4 / 1.2 - 1, 'relative_step_error': 1 / 11, 'margins': margins(f), 'poles': [[x.real, x.imag] for x in np.roots(den + np.array([0, 0, 0, 10]))]}
    w = np.geomspace(0.001, 100, 5000)
    z = f(w)
    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    ax.plot(z.real, z.imag, label='omega > 0')
    ax.plot(z.real, -z.imag, ls='--', label='omega < 0')
    ax.plot(-1, 0, 'rx', ms=8)
    ax.set(xlim=(-1.3, 10.5), ylim=(-7, 7), xlabel='Re L(j omega)', ylabel='Im L(j omega)')
    ax.axhline(0, color='k', lw=0.6)
    ax.axvline(0, color='k', lw=0.6)
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / 'P9_5_nyquist.png')
    plt.close(fig)
    np.savetxt(HERE / 'P9_5_nyquist.csv', np.c_[w, z.real, z.imag], delimiter=',', header='omega,real,imag', comments='')
    base = tf([2], np.polymul([0.2, 0.1], np.polymul([1, 0], [0.4, 1])))
    wc = optimize.brentq(lambda x: np.arctan(2 * x) + np.arctan(0.4 * x) - np.pi / 3, 0.001, 10)
    k = 1 / abs(base(wc))
    f = lambda w: k * base(w)
    ze = np.sin(np.pi / 6) / (2 * np.sqrt(np.cos(np.pi / 6)))
    wn = wc / np.sqrt(np.sqrt(1 + 4 * ze ** 4) - 2 * ze ** 2)
    n = [2 * k]
    d = np.polymul([0.2, 0.1], np.polymul([1, 0], [0.4, 1]))
    d[-1] += 2 * k
    RESULT['P9.10'] = {'signed_K1': -k, 'gain_magnitude': k, 'wc': wc, 'second_order_zeta': ze, 'second_order_wn': wn, 'estimate_Mr': 1 / (2 * ze * np.sqrt(1 - ze ** 2)), 'estimate_wr': wn * np.sqrt(1 - 2 * ze ** 2), 'estimate_bandwidth': wn * np.sqrt(1 - 2 * ze ** 2 + np.sqrt(2 - 4 * ze ** 2 + 4 * ze ** 4)), 'estimate_overshoot_pct': 100 * np.exp(-np.pi * ze / np.sqrt(1 - ze ** 2)), 'estimate_settling': 4 / (ze * wn), 'exact_frequency': closed_metrics(f), 'exact_step': step('P9_10_step', n, d, np.linspace(0, 80, 40001)), 'margins': margins(f)}
    bode('P9_10_bode', [f], ['L'], np.geomspace(0.01, 100, 4000), [wc])
    base = lambda w: np.exp(-0.2j * np.asarray(w)) / (1j * np.asarray(w) * (1 + 0.1j * np.asarray(w)))
    kcritw = optimize.brentq(lambda x: np.arctan(0.1 * x) + 0.2 * x - np.pi / 2, 0.01, 20)
    kcrit = 1 / abs(base(kcritw))
    target = 10 ** 0.1
    K = optimize.brentq(lambda x: closed_metrics(lambda w: x * base(w))['Mr'] - target, 0.05, kcrit * 0.999)
    f = lambda w: K * base(w)
    met = closed_metrics(f)
    mar = margins(f, hi=1000)
    zeMr = np.sqrt((1 - np.sqrt(1 - 1 / met['Mr'] ** 2)) / 2)
    pm = mar['gain_crossovers'][0]['PM_deg'] * np.pi / 180
    zePM = np.sin(pm) / (2 * np.sqrt(np.cos(pm)))
    RESULT['P9.14'] = {'K': K, 'critical_K': kcrit, 'critical_frequency': kcritw, **met, 'PM_deg': pm * 180 / np.pi, 'zeta_from_Mr': zeMr, 'zeta_from_PM': zePM}
    theta = np.linspace(-359, -1, 600)
    db = np.linspace(-35, 25, 500)
    TH, DB = np.meshgrid(theta, db)
    Z = 10 ** (DB / 20) * np.exp(1j * np.deg2rad(TH))
    M = 20 * np.log10(abs(Z / (1 + Z)))
    fig, ax = plt.subplots(figsize=(7.1, 4.4))
    cs = ax.contour(TH, DB, M, levels=[-3, 0, 1, 2, 3, 6], colors='grey', linewidths=0.7)
    ax.clabel(cs, levels=[-3, 1, 2, 3, 6], fmt='%g dB', fontsize=8)
    w = np.geomspace(0.02, 8, 4000)
    zz = f(w)
    ph = (-np.pi / 2 - np.arctan(0.1 * w) - 0.2 * w) * 180 / np.pi
    ax.plot(ph, 20 * np.log10(abs(zz)), lw=2, label=f'K={K:.5f}')
    ax.plot(-180, 0, 'rx')
    ax.set(xlim=(-280, -80), ylim=(-30, 25), xlabel='Open-loop phase (deg)', ylabel='Open-loop magnitude (dB)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / 'P9_14_nichols.png')
    plt.close(fig)
    np.savetxt(HERE / 'P9_14_nichols.csv', np.c_[w, ph, 20 * np.log10(abs(zz)), abs(zz / (1 + zz))], delimiter=',', header='omega,phase_deg,magnitude_dB,closed_magnitude', comments='')
    K = np.sqrt(2)
    ze = 1 / (2 * np.sqrt(K))
    phi = np.pi / 4
    zepm = np.sin(phi) / (2 * np.sqrt(np.cos(phi)))
    RESULT['P9.28'] = {'Kp': K, 'wc': 1, 'zeta': ze, 'zeta_from_PM': zepm, 'overshoot_PM_pct': 100 * np.exp(-np.pi * zepm / np.sqrt(1 - zepm * zepm)), **step('P9_28_step', [K], [1, 1, K], np.linspace(0, 16, 12001))}
    fs = []
    for wn2 in [15267, 9500]:
        n = np.polymul([775, 7750], [305.3])
        n = np.polymul(n, [1, 5])
        d = np.polymul(np.polymul([1, 0], [1, 2]), np.polymul([1, 10, wn2], [1, 1]))
        f = tf(n, d)
        fs.append(f)
        D = np.polyadd(d, n)
        RESULT['AP9.1_' + str(wn2)] = {'margins': margins(f), 'poles': [[x.real, x.imag] for x in np.roots(D)]}
    bode('AP9_1_bode', fs, ['wn^2=15267', 'wn^2=9500'], np.geomspace(0.03, 1000, 10000))
    (HERE / 'numeric_results.json').write_text(json.dumps(RESULT, ensure_ascii=False, indent=2), encoding='utf8')
    print(json.dumps(RESULT, ensure_ascii=True, indent=2))
if __name__ == '__main__':
    run()
