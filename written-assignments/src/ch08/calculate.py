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
import os, json
BASE = DATA
import numpy as np
from scipy import optimize, signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out = {}

def solvelevel(fun, lo, hi, level=1):
    return optimize.brentq(lambda x: abs(fun(1j * x)) - level, lo, hi, xtol=1e-11)

def bode(name, fun, w, title, asym=None, marks=None):
    g = fun(1j * w)
    mag = 20 * np.log10(abs(g))
    phase = np.degrees(np.unwrap(np.angle(g)))
    np.savetxt(BASE / f'{name}.csv', np.c_[w, g.real, g.imag, mag, phase], delimiter=',', header='omega_rad_s,re,im,mag_db,phase_deg', comments='')
    fig, ax = plt.subplots(2, 1, figsize=(7, 5.3), sharex=True, layout='constrained')
    ax[0].semilogx(w, mag, label='Exact')
    ax[1].semilogx(w, phase, label='Exact')
    if asym is not None:
        ax[0].semilogx(w, asym[0], ls='--', label='Asymptote')
        ax[1].semilogx(w, asym[1], ls='--', label='Phase approximation')
    if marks:
        for wm, label in marks:
            gm = fun(1j * wm)
            ax[0].plot(wm, 20 * np.log10(abs(gm)), 'o', ms=4)
            ax[1].plot(wm, np.degrees(np.angle(gm)), 'o', ms=4, label=label)
            for a in ax:
                a.axvline(wm, color='gray', lw=0.7, ls=':')
    ax[0].axhline(0, color='gray', lw=0.6)
    ax[0].set(ylabel='Magnitude (dB)', title=title)
    ax[1].set(xlabel='Frequency (rad/s)', ylabel='Phase (deg)')
    for a in ax:
        a.grid(True, which='both', alpha=0.25)
        a.legend(fontsize=8)
    fig.savefig(FIGURES / f'{name}.png', dpi=200)
    plt.close(fig)
F = lambda s: 20 / (s * (1 + s / 8) * (1 + s / 120))
wc = solvelevel(F, 0.1, 100)
out['E8.6'] = {'wc': wc, 'phase_wc': float(np.angle(F(1j * wc), deg=True)), 'magnitude_residual': float(abs(F(1j * wc)) - 1)}
bode('E8_6', F, np.logspace(-1, 4, 3000), 'E8.6: arm-control Bode diagram', marks=[(wc, f'wc = {wc:.4f} rad/s')])
F2 = lambda s: 20 / (s * (s + 1.4) * (s + 10))
wc2 = solvelevel(F2, 0.1, 10)
out['E8.14'] = {'wc': wc2, 'phase_wc': float(np.angle(F2(1j * wc2), deg=True)), 'magnitude_residual': float(abs(F2(1j * wc2)) - 1)}
bode('E8_14', F2, np.logspace(-2, 3, 3000), 'E8.14: non-unit feedback loop Bode diagram', marks=[(wc2, f'wc = {wc2:.4f} rad/s')])
polars = [lambda s: 1 / ((1 + 0.25 * s) * (1 + 3 * s)), lambda s: 5 * (s * s + 1.4 * s + 1) / (s - 1) ** 2, lambda s: (s - 8) / (s * s + 6 * s + 8), lambda s: 20 * (s + 8) / (s * (s + 2) * (s + 4))]
limits = [(-0.25, 1.1, -0.75, 0.15), (-5, 6, -5, 5), (-1.1, 0.5, -0.25, 0.9), (-14, 1, -10, 4)]
fig, axes = plt.subplots(2, 2, figsize=(10, 8), layout='constrained')
out['P8.1'] = {}
for i, (fun, lim) in enumerate(zip(polars, limits)):
    w = np.logspace(-5, 5, 16000)
    g = fun(1j * w)
    ax = axes.flat[i]
    ax.plot(g.real, g.imag, label='0 < omega < infinity')
    ax.axhline(0, color='gray', lw=0.6)
    ax.axvline(0, color='gray', lw=0.6)
    ax.grid(True, alpha=0.3)
    ax.set(xlabel='Re L(j omega)', ylabel='Im L(j omega)', title=f'P8.1({chr(97 + i)})', xlim=lim[:2], ylim=lim[2:])
    for om in [0.2, 1, 5] if i != 3 else [1.5, 5, 20]:
        a = fun(1j * om)
        b = fun(1j * om * 1.12)
        ax.annotate('', xy=(b.real, b.imag), xytext=(a.real, a.imag), arrowprops={'arrowstyle': '->', 'color': 'C0'})
    if i == 3:
        ax.axvline(-12.5, color='gray', ls='--', lw=0.8, label='omega -> 0: Re = -12.5')
        zz = ax.inset_axes([0.5, 0.42, 0.44, 0.3])
        zz.plot(g.real, g.imag)
        zz.axhline(0, color='gray', lw=0.5)
        zz.axvline(0, color='gray', lw=0.5)
        zz.set(xlim=(-1.6, 0.05), ylim=(-0.12, 0.12), title='Near-origin detail')
        zz.tick_params(labelsize=7)
        zz.title.set_fontsize(8)
        zz.grid(True, alpha=0.2)
        a = fun(1j * np.sqrt(32))
        zz.plot(a.real, a.imag, 'o', ms=3)
    ax.legend(fontsize=7)
    np.savetxt(BASE / f'P8_1{chr(97 + i)}.csv', np.c_[w, g.real, g.imag], delimiter=',', header='omega_rad_s,re,im', comments='')
fig.savefig(FIGURES / 'P8_1.png', dpi=200)
plt.close(fig)
w = np.logspace(-3, 3, 4000)
forms = [((1 - 0.75 * w * w) / ((1 - 0.75 * w * w) ** 2 + (3.25 * w) ** 2), -3.25 * w / ((1 - 0.75 * w * w) ** 2 + (3.25 * w) ** 2)), (5 * ((1 - w * w) ** 2 - 2.8 * w * w) / (1 + w * w) ** 2, 17 * w * (1 - w * w) / (1 + w * w) ** 2), ((14 * w * w - 64) / ((8 - w * w) ** 2 + 36 * w * w), w * (56 - w * w) / ((8 - w * w) ** 2 + 36 * w * w)), (-20 * (w * w + 40) / ((w * w + 4) * (w * w + 16)), 40 * (w * w - 32) / (w * (w * w + 4) * (w * w + 16)))]
out['P8.1']['expression_max_errors'] = [float(max(abs(re + 1j * im - fun(1j * w)))) for fun, (re, im) in zip(polars, forms)]
theta = lambda w, c: np.clip(45 * np.log10(10 * w / c), 0, 90)
w = np.logspace(-2, 4, 4000)
K = 10 ** 0.6
fa = lambda s: K * (1 + s) / (s * (1 + s / 10) ** 2)
fb = lambda s: s / ((1 + s / 10) * (1 + s / 50))
ma = 20 * np.log10(K) - 20 * np.log10(w) + 20 * np.maximum(0, np.log10(w)) - 40 * np.maximum(0, np.log10(w / 10))
pa = -90 + theta(w, 1) - 2 * theta(w, 10)
mb = 20 * np.log10(w) - 20 * np.maximum(0, np.log10(w / 10)) - 20 * np.maximum(0, np.log10(w / 50))
pb = 90 - theta(w, 10) - theta(w, 50)
fig, axs = plt.subplots(2, 2, figsize=(10, 6.5), sharex='col', layout='constrained')
for j, (fun, mag, ph, label) in enumerate([(fa, ma, pa, 'P8.6(a): repeated real-pole model'), (fb, mb, pb, 'P8.6(b)')]):
    g = fun(1j * w)
    pp = np.unwrap(np.angle(g)) * 180 / np.pi
    axs[0, j].semilogx(w, 20 * np.log10(abs(g)), label='Exact')
    axs[0, j].semilogx(w, mag, '--', label='Asymptote')
    axs[0, j].set(title=label, ylabel='Magnitude (dB)')
    axs[1, j].semilogx(w, pp, label='Exact')
    axs[1, j].semilogx(w, ph, '--', label='Phase approximation')
    axs[1, j].set(xlabel='Frequency (rad/s)', ylabel='Phase (deg)')
    np.savetxt(BASE / f'P8_6{chr(97 + j)}.csv', np.c_[w, 20 * np.log10(abs(g)), mag, pp, ph], delimiter=',', header='omega_rad_s,mag_db,asymptote_db,phase_deg,phase_approx_deg', comments='')
for ax in axs.flat:
    ax.grid(True, which='both', alpha=0.25)
    ax.legend(fontsize=7)
fig.savefig(FIGURES / 'P8_6.png', dpi=200)
plt.close(fig)
out['P8.6'] = {'K_a': K, 'b_corners': [10, 50]}
Fg = lambda s: 1000000000.0 * s / ((s + 1000) * (s + 10000000.0))
peakw = 100000.0
peak = abs(Fg(1j * peakw))
abslevel = 10 ** (-3 / 20)
out['P8.14'] = {'midband_gain': 100, 'peak': peak, 'peak_frequency': peakw, 'absolute_minus3db': [solvelevel(Fg, 1, 1000, abslevel), solvelevel(Fg, 10000000.0, 100000000000.0, abslevel)], 'relative_peak_halfpower': [solvelevel(Fg, 100, 100000.0, peak / np.sqrt(2)), solvelevel(Fg, 100000.0, 100000000.0, peak / np.sqrt(2))]}
bode('P8_14', Fg, np.logspace(0, 11, 5000), 'P8.14: band-pass amplifier', marks=[(1000, 'Low corner = 1000 rad/s'), (10000000.0, 'High corner = 1e7 rad/s')])
mr = 10 ** (12 / 20)
z = np.sqrt((1 - np.sqrt(1 - 1 / mr ** 2)) / 2)
wn = 0.9
wb = wn * np.sqrt(1 - 2 * z * z + np.sqrt(2 - 4 * z * z + 4 * z ** 4))
mp = 100 * np.exp(-np.pi * z / np.sqrt(1 - z * z))
ts = 4 / (z * wn)
out['P8.24'] = {'read_peak_db': 12, 'read_wn': wn, 'Mr': mr, 'zeta': z, 'denominator': [1, 2 * z * wn, wn * wn], 'bandwidth_rad_s': wb, 'overshoot_percent': mp, 'ts_2percent_estimate_s': ts}
tt = np.linspace(0, 65, 8000)
t, y = signal.step(signal.TransferFunction([wn * wn], [1, 2 * z * wn, wn * wn]), T=tt)
np.savetxt(BASE / 'P8_24_step.csv', np.c_[t, y], delimiter=',', header='time_s,y', comments='')
fig, ax = plt.subplots(figsize=(7, 3.8), layout='constrained')
ax.plot(t, y, label=f'wn = 0.9 rad/s, zeta = {z:.4f}')
ax.axhline(1, color='gray', ls='--')
ax.fill_between(t, 0.98, 1.02, color='gray', alpha=0.2, label='2% band')
ax.set(xlabel='Time (s)', ylabel='Output y(t)', title='P8.24: step response of the second-order approximation')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)
fig.savefig(FIGURES / 'P8_24.png', dpi=200)
plt.close(fig)
(BASE / 'numeric_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding='utf-8')
print(json.dumps(out, indent=2))
