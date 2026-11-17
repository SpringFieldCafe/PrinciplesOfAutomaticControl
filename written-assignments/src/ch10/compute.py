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
import os, json
H = DATA
import numpy as np
from scipy import signal, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.25, 'savefig.dpi': 210})
R = {}

def tf(n, d):
    return lambda w: np.polyval(n, 1j * np.asarray(w)) / np.polyval(d, 1j * np.asarray(w))

def cross(fun, lo=1e-05, hi=100000.0):
    w = np.geomspace(lo, hi, 30000)
    v = fun(w)
    ix = np.where(v[:-1] * v[1:] < 0)[0]
    return [float(optimize.brentq(fun, w[i], w[i + 1], xtol=1e-12)) for i in ix]

def margin(f):
    w = cross(lambda x: abs(f(x)) - 1)
    return [{'wc': x, 'PM': float(180 + np.angle(f(x), deg=True))} for x in w]

def metrics(t, y, final=1):
    q = y / final
    ids = np.where(abs(q - 1) > 0.02)[0]
    settle = float(t[ids[-1] + 1]) if len(ids) and ids[-1] + 1 < len(t) else None

    def timeat(v):
        ix = np.where(q >= v)[0]
        if not len(ix):
            return None
        i = ix[0]
        return float(t[i]) if i == 0 else float(t[i - 1] + (v - q[i - 1]) * (t[i] - t[i - 1]) / (q[i] - q[i - 1]))
    t10 = timeat(0.1)
    t90 = timeat(0.9)
    return {'overshoot_pct': float(max(0, (max(q) - 1) * 100)), 'settling_2pct': settle, 'peak_time': float(t[np.argmax(q)]), 'rise_10_90': t90 - t10 if t90 is not None else None, 'rise_0_100': timeat(1)}

def response(n, d, end=20, npts=20001):
    t, y = signal.step((n, d), T=np.linspace(0, end, npts))
    return (t, y)

def stepfig(name, items, end):
    fig, ax = plt.subplots(figsize=(7.1, 4))
    cols = []
    data = {}
    for label, n, d in items:
        t, y = response(n, d, end)
        ax.plot(t, y, label=label)
        cols.append(y)
        data[label] = metrics(t, y)
    ax.axhline(1, color='k', lw=0.7)
    ax.axhspan(0.98, 1.02, alpha=0.1, color='grey')
    ax.set(xlabel='Time (s)', ylabel='Output')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / (name + '.png'))
    plt.close(fig)
    np.savetxt(H / (name + '.csv'), np.column_stack([t, *cols]), delimiter=',', header='time_s,' + ','.join((x[0] for x in items)), comments='')
    return data

def freqmetrics(f):
    w = np.geomspace(1e-05, 10000.0, 40000)
    a = abs(f(w) / (1 + f(w)))
    i = np.argmax(a)
    fit = optimize.minimize_scalar(lambda x: -abs(f(np.exp(x)) / (1 + f(np.exp(x)))), bounds=(np.log(w[max(0, i - 2)]), np.log(w[min(len(w) - 1, i + 2)])), method='bounded')
    b = cross(lambda x: abs(f(x) / (1 + f(x))) - 1 / np.sqrt(2))[0]
    return {'Mr': float(-fit.fun), 'wr': float(np.exp(fit.x)), 'bandwidth': b}

def bodefig(name, fs, labels):
    w = np.geomspace(0.01, 1000, 6000)
    fig, ax = plt.subplots(2, 1, figsize=(7.1, 5.4), sharex=True)
    cols = []
    for fun, label in zip(fs, labels):
        z = fun(w)
        a = 20 * np.log10(abs(z))
        ph = np.unwrap(np.angle(z)) * 180 / np.pi
        cols.extend([a, ph])
        ax[0].semilogx(w, a, label=label)
        ax[1].semilogx(w, ph)
    ax[0].axhline(0, color='k', lw=0.6)
    ax[1].axhline(-180, color='k', lw=0.6)
    ax[0].set_ylabel('Magnitude (dB)')
    ax[1].set_ylabel('Phase (deg)')
    ax[1].set_xlabel('Frequency (rad/s)')
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(FIGURES / (name + '.png'))
    plt.close(fig)
    np.savetxt(H / (name + '.csv'), np.column_stack([w, *cols]), delimiter=',', header='omega,' + ','.join((f'{label}_{v}' for label in labels for v in ['dB', 'deg'])), comments='')

def poles(d):
    return [[float(p.real), float(p.imag)] for p in np.roots(d)]

def delay_response(K, b=None, end=40, h=0.002):
    n = [K] if b is None else [K, 2 * K]
    d = [1, 4, 3] if b is None else np.polymul([1, 4, 3], [1, b])
    A, B, C, D = signal.tf2ss(n, d)
    B = B[:, 0]
    C = C[0]
    nt = round(end / h)
    delay = round(0.5 / h)
    y = np.zeros(nt + 1)
    x = np.zeros(A.shape[0])
    for j in range(nt):

        def inp(frac):
            at = j + frac - delay
            if at < 0:
                return 0.0
            ix = int(np.floor(at))
            rem = at - ix
            old = (1 - rem) * y[ix] + rem * y[ix + 1]
            return 1 - old
        k1 = A @ x + B * inp(0)
        k2 = A @ (x + h * k1 / 2) + B * inp(0.5)
        k3 = A @ (x + h * k2 / 2) + B * inp(0.5)
        k4 = A @ (x + h * k3) + B * inp(1)
        x = x + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        y[j + 1] = C @ x
    t = np.arange(nt + 1) * h
    fin = K / (3 + K) if b is None else 2 * K / (3 * b + 2 * K)
    return (t, y, metrics(t, y, fin), fin)

def main():
    K = np.log(50)
    R['E10.1'] = {'a': 2.0, 'K': K, 'exact_Ts': np.log(50) / K, 'step': stepfig('E10_1_step', [('PI design', [K], [1, K])], 2.5)}
    n = [10, 40]
    d = [1, 15.2, 153, 40, 40]
    pp = np.roots(d)
    dom = max(pp, key=lambda p: p.real)
    wn = abs(dom)
    ze = -dom.real / wn
    R['E10.6'] = {'closed_numerator': n, 'closed_denominator': d, 'poles': poles(d), 'dominant_zeta': ze, 'dominant_wn': wn, 'estimate_overshoot_pct': 100 * np.exp(-np.pi * ze / np.sqrt(1 - ze * ze)), 'estimate_settling': 4 / (ze * wn), 'step': stepfig('E10_6_step', [('Full fourth-order model', n, d)], 65)}

    def piPO(kp):
        sig = (1 + kp) / 2
        wd = np.sqrt(5 - sig * sig)
        err = lambda t: np.exp(-sig * t) * (np.cos(wd * t) + (1 - sig) * np.sin(wd * t) / wd)
        best = optimize.minimize_scalar(err, bounds=(0.01, 3), method='bounded')
        return -100 * best.fun
    kp = optimize.brentq(lambda k: piPO(k) - 5, 2.5, 3.3)
    R['E10.9'] = {'Kp': kp, 'Ki': 5.0, 'Kv': 5.0, 'poles': poles([1, 1 + kp, 5]), 'step': stepfig('E10_9_step', [('PI design', [kp, 5], [1, 1 + kp, 5])], 8)}
    designs = []
    items = []
    fs = []
    for label, r in [('Bode design', 20.0), ('Root-locus design', 12.0)]:
        wn = 4
        zeta = 0.6
        p = r + 2 * zeta * wn
        A = wn * wn + 2 * zeta * wn * r
        z = r * wn * wn / A
        n = [A, A * z]
        d = [1, p, A, A * z]
        f = tf(n, [1, p, 0, 0])
        fs.append(f)
        items.append((label, n, d))
        designs.append({'method': label, 'p': p, 'A': A, 'z': z, 'poles': poles(d), 'margin': margin(f), 'frequency': freqmetrics(f)})
    R['P10.1'] = {'designs': designs, 'step': stepfig('P10_1_step', items, 4)}
    bodefig('P10_1_bode', fs, ['Bode design', 'Root-locus design'])
    vals = np.geomspace(0.1, 5000, 1400)
    rr = np.array([np.roots([1, 16.8, a, a * (192 / 73.6)]) for a in vals])
    fig, ax = plt.subplots(figsize=(6.8, 4))
    for i in range(3):
        ax.scatter(rr[:, i].real, rr[:, i].imag, s=1, color='C0')
    ax.plot([-2.4, -2.4], [3.2, -3.2], 'rx', ms=8, label='Selected closed-loop poles')
    ax.plot(-12, 0, 'rx', ms=8)
    ax.set(xlim=(-20, 2), ylim=(-12, 12), xlabel='Real part (1/s)', ylabel='Imaginary part (rad/s)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / 'P10_1_rootlocus.png')
    plt.close(fig)
    np.savetxt(H / 'P10_1_rootlocus.csv', np.column_stack([vals, rr.real, rr.imag]), delimiter=',', header='A,real1,real2,real3,imag1,imag2,imag3', comments='')
    ze = -np.log(0.1) / np.sqrt(np.pi ** 2 + np.log(0.1) ** 2)
    wn = 1 / 0.3 / ze
    kp = wn * wn / 25
    ki = 0.15 * kp
    d = [1, 20 / 3, wn * wn]
    R['P10.5'] = {'Ka': 1.0, 'Kp': kp, 'Ki': ki, 'zeta': ze, 'wn': wn, 'estimate_settling': 1.2, 'step': stepfig('P10_5_step', [('PI design', [wn * wn], d)], 3)}
    base = tf([4], [1, 3, 2, 0])
    wc = optimize.brentq(lambda w: np.arctan(w) + np.arctan(w / 2) - np.deg2rad(40), 0.01, 2)
    K = 1 / abs(base(wc))
    fa = lambda w: K * base(w)
    z = 0.05

    def lagres(v):
        pole, w = np.exp(v)
        g = (1 + 1j * w / z) / (1 + 1j * w / pole)
        val = g * base(w)
        return [np.log(abs(val)), np.angle(val) + np.deg2rad(130)]
    lr = optimize.root(lagres, np.log([0.012, 0.45]))
    pole, w = np.exp(lr.x)
    n = np.polymul([4], [pole / z, pole])
    d = np.polymul([1, 3, 2, 0], [1, pole])
    fc = tf(n, d)
    R['P10.13'] = {'part_a': {'K': K, 'wc': wc, 'margin': margin(fa), **freqmetrics(fa)}, 'part_b': {'lag_zero': z, 'lag_pole': pole, 'Kv': 2.0, 'margin': margin(fc), 'poles': poles(np.polyadd(d, n)), **freqmetrics(fc), 'residual': lagres(lr.x)}}
    bodefig('P10_13_bode', [fa, fc], ['Proportional', 'Lag compensated'])

    def poK(k):
        return delay_response(k, end=22, h=0.004)[2]['overshoot_pct']
    maxK = optimize.brentq(lambda k: poK(k) - 30, 1, 10, xtol=1e-07)
    K = 1.0
    ta, ya, ma, fina = delay_response(K, end=40, h=0.002)
    trials = []
    for b in [0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]:
        t, y, m, fin = delay_response(K, b, end=45, h=0.004)
        trials.append({'b': b, 'steady_error': 1 - fin, **m})
    b = 1 / 15
    tb, yb, mb, finb = delay_response(K, b, end=45, h=0.002)
    converge = []
    for hh in [0.004, 0.002, 0.001]:
        _, _, m2, ff = delay_response(K, b, end=45, h=hh)
        converge.append({'h': hh, 'final': ff, **m2})
    R['P10.19'] = {'part_a': {'K': K, 'final': fina, 'steady_error': 1 - fina, **ma}, 'part_b': {'K': K, 'b': b, 'final': finb, 'steady_error': 1 - finb, **mb}, 'lag_trials': trials, 'step_convergence': converge, 'uncompensated_K_for_30pct': maxK}
    fig, ax = plt.subplots(figsize=(7.1, 4))
    ax.plot(ta, ya, label=f'Uncompensated K={K:.4f}')
    ax.plot(tb, yb, label=f'Lag b={b}')
    ax.axhline(1, color='k', lw=0.7)
    ax.set(xlabel='Time (s)', ylabel='Output')
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURES / 'P10_19_step.png')
    plt.close(fig)
    np.savetxt(H / 'P10_19_step.csv', np.c_[tb, np.interp(tb, ta, ya), yb], delimiter=',', header='time_s,uncompensated,compensated', comments='')
    base = tf([20], np.polymul([1, 0], np.polymul([0.5, 1], [1 / 6, 1])))
    wc = 10.0
    delta = (np.deg2rad(-135) - np.angle(base(wc))) / 2
    if delta < 0:
        delta += np.pi
    M = np.sqrt(1 / abs(base(wc)))
    y = (np.cos(delta) - 1 / M) / np.sin(delta)
    x = M * (np.sin(delta) + y * np.cos(delta))
    z = wc / x
    p = wc / y
    n = np.polymul([20], np.polymul([1 / z, 1], [1 / z, 1]))
    d = np.polymul(np.polymul([1, 0], np.polymul([0.5, 1], [1 / 6, 1])), np.polymul([1 / p, 1], [1 / p, 1]))
    fc = tf(n, d)
    R['P10.20'] = {'K': 20.0, 'lead_zero': z, 'lead_pole': p, 'phase_each_deg': np.degrees(delta), 'margin': margin(fc), 'Kv': 20.0, **freqmetrics(fc), 'poles': poles(np.polyadd(d, n))}
    bodefig('P10_20_bode', [base, fc], ['Uncompensated', 'Two identical leads'])
    wc = 3.0

    def lagleadres(v):
        a, b = np.exp(v)
        s = 1j * wc
        g = (1 + s / (10 * a)) * (1 + s / b) / ((1 + s / a) * (1 + s / (10 * b)))
        val = base(wc) * g
        return [np.log(abs(val)), np.angle(val) + np.deg2rad(135)]
    fit = optimize.root(lagleadres, np.log([0.1, 1.0]))
    a, b = np.exp(fit.x)
    n = np.polymul([20], np.polymul([1 / (10 * a), 1], [1 / b, 1]))
    d = np.polymul(np.polymul([1, 0], np.polymul([0.5, 1], [1 / 6, 1])), np.polymul([1 / a, 1], [1 / (10 * b), 1]))
    fc = tf(n, d)
    R['P10.22'] = {'K': 20.0, 'a': a, 'b': b, 'Kv': 20.0, 'residual': lagleadres(fit.x), 'margin': margin(fc), **freqmetrics(fc), 'poles': poles(np.polyadd(d, n))}
    bodefig('P10_22_bode', [base, fc], ['Uncompensated', 'Lag-lead compensated'])
    kl = 35 - 4 * np.sqrt(76)
    ku = 35 + 4 * np.sqrt(76)

    def modalPO(k):
        n = [k, 2 * k, 20 * k]
        d = [1, 4, 5 + k, 2 + 2 * k, 20 * k]
        rr, pp, _ = signal.residue(n, np.polymul(d, [1, 0]))
        ix = abs(pp) > 1e-08
        rr = rr[ix]
        pp = pp[ix]
        maxT = min(200000, 12 / -max(pp.real))
        tt = np.unique(np.r_[np.linspace(0, min(80, maxT), 18000), np.geomspace(1e-05, maxT, 10000)])
        yy = 1 + np.real(np.exp(np.outer(tt, pp)) @ rr)
        i = np.argmax(yy)
        lo = tt[max(i - 1, 0)]
        hi = tt[min(i + 1, len(tt) - 1)]
        fit = optimize.minimize_scalar(lambda t: -(1 + np.real(np.exp(t * pp) @ rr)), bounds=(lo, hi), method='bounded')
        return max(0, (-fit.fun - 1) * 100)
    k1 = np.r_[np.geomspace(0.001, kl * 0.98, 70), np.linspace(kl * 0.98, kl * (1 - 1e-05), 35)]
    k2 = np.linspace(ku * (1 + 1e-05), 100, 120)
    po1 = [modalPO(k) for k in k1]
    po2 = [modalPO(k) for k in k2]
    fig, ax = plt.subplots(figsize=(7.2, 4.1))
    ax.plot(k1, po1, 'C0')
    ax.plot(k2, po2, 'C0')
    ax.axvspan(kl, ku, color='grey', alpha=0.18, label='Unstable: no finite steady-state overshoot')
    ax.set(xlabel='Gain K', ylabel='Overshoot (%)', xlim=(0, 100))
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURES / 'P10_43_overshoot.png')
    plt.close(fig)
    np.savetxt(H / 'P10_43_overshoot.csv', np.c_[np.r_[k1, k2], np.r_[po1, po2]], delimiter=',', header='K,overshoot_pct', comments='')
    R['P10.43'] = {'stable_lower_upper': kl, 'stable_high_lower': ku, 'sample_poles': {str(k): poles([1, 4, 5 + k, 2 + 2 * k, 20 * k]) for k in [0.1, 1, 10, 40, 80, 100]}, 'endpoint_overshoots': [po1[0], po1[-1], po2[0], po2[-1]]}
    K1 = 6.4
    K2 = 0.04375
    R['AP10.4'] = {'K1': K1, 'K2': K2, 'wn': 8.0, 'zeta': 0.8, 'estimate_settling': 0.625, 'step': stepfig('AP10_4_step', [('Speed-feedback design', [64], [1, 12.8, 64])], 1.5)}
    xs = np.roots([1, -38, 205, -444, 336])
    choices = []
    for xx in xs:
        if abs(xx.imag) < 1e-08 and xx.real > 0:
            x = xx.real
            K = -76 * x * x + 592 * x - 448
            d = [1, 38, 296, 448 + K, 3 * K]
            if K > 0 and max(np.roots(d).real) < 0:
                choices.append((K, x, d))
    K, x, d = min(choices, key=lambda a: a[0])
    n = [K, 3 * K]

    def filteredPO(p):
        tn, yn = response(np.polymul(n, [p]), np.polymul(d, [1, p]), 8, 18001)
        return metrics(tn, yn)['overshoot_pct']
    p = optimize.brentq(lambda p: filteredPO(p) - 1, 0.1, 3, xtol=1e-09)
    ns = np.polymul(n, [p])
    ds = np.polymul(d, [1, p])
    r3n = np.polymul(n, [3])
    r3d = np.polymul(d, [1, 3])
    R['AP10.7'] = {'K': K, 'dominant_sigma': x, 'all_choices': [(a, b) for a, b, _ in choices], 'poles': poles(d), 'p_for_1pct': p, 'step': stepfig('AP10_7_step', [('No prefilter', n, d), ('Prefilter p=3', r3n, r3d), (f'Prefilter p={p:.4f}', ns, ds)], 7)}
    (H / 'numeric_results.json').write_text(json.dumps(R, ensure_ascii=False, indent=2), encoding='utf8')
    print(json.dumps(R, ensure_ascii=True, indent=2))
if __name__ == '__main__':
    main()
