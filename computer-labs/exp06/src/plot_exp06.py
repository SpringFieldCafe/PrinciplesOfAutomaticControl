"""Plot the actual CSV results produced by run_exp06 (NumPy/Matplotlib).

Run with the existing Anaconda Python. All outputs and caches stay in the
project's ignored .local/exp06 directory; no numerical fitting is done here.
"""
from pathlib import Path
import os

ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.local' / 'exp06'
OUT = PRIVATE / 'results'
for folder in (OUT, PRIVATE / 'mplconfig', PRIVATE / 'tmp'):
    folder.mkdir(parents=True, exist_ok=True)
os.environ['MPLCONFIGDIR'] = str(PRIVATE / 'mplconfig')
os.environ['TMP'] = os.environ['TEMP'] = str(PRIVATE / 'tmp')

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, MultipleLocator

plt.rcParams.update({
    'font.family': ['Microsoft YaHei', 'DejaVu Sans'], 'font.size': 10,
    'axes.unicode_minus': False, 'axes.spines.top': False,
    'axes.spines.right': False, 'axes.titlepad': 12,
    'svg.fonttype': 'none', 'savefig.facecolor': 'white',
})
BLUE, ORANGE, GREEN, GRAY = '#17649b', '#ca6b22', '#20856b', '#6b7280'
FORMULA = r'$G(s)=10/[s(s+2)(s+5)]$'


def load(name):
    data = np.genfromtxt(OUT / (name + '.csv'), delimiter=',', names=True,
                         encoding='utf-8-sig')
    if data.ndim != 1 or len(data) < 2:
        raise ValueError(f'{name}: expected at least two CSV rows')
    for col in data.dtype.names:
        if not np.all(np.isfinite(data[col])):
            raise ValueError(f'{name}: nonfinite values in {col}')
    if np.any(data['w_rad_s'] <= 0) or np.any(np.diff(data['w_rad_s']) <= 0):
        raise ValueError(f'{name}: frequency must be positive and increasing')
    return data


def decorate(ax, ylabel, xlabel=True):
    ax.set_xscale('log')
    ax.set_ylabel(ylabel)
    if xlabel:
        ax.set_xlabel(r'角频率 $\omega$ (rad/s)')
    ax.grid(which='major', color='#bec6ce', lw=.65, alpha=.65)
    ax.grid(which='minor', color='#d8dfe5', lw=.4, alpha=.5)
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.xaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax.set_xlim(.1, 100)


def pair(title, note=None, height=7.8):
    fig, ax = plt.subplots(2, 1, figsize=(10.6, height), sharex=True)
    fig.subplots_adjust(left=.1, right=.96, top=.87, bottom=.12, hspace=.19)
    fig.suptitle(title, fontsize=15, y=.965)
    fig.text(.5, .915, FORMULA, ha='center', fontsize=12, color=GRAY)
    decorate(ax[0], '幅值 (dB)', False)
    decorate(ax[1], '相位 (°)')
    if note:
        fig.text(.1, .035, note, fontsize=9, color=GRAY)
    return fig, ax


def save(fig, name):
    for ext in ('png', 'svg'):
        fig.savefig(OUT / f'{name}.{ext}', dpi=240)
    plt.close(fig)
    print(name)


def asymptotes(a, sketch=False):
    note = ('程序绘制的笔绘风格：坐标、刻度及关键点固定；曲线采用轻微视觉笔触。'
            if sketch else '标准绘图：幅值按分段斜率构造；相位按各一阶因子的两十倍频程近似相加。')
    fig, ax = pair('任务 A｜Bode 渐近线' + (' · 笔绘风格' if sketch else ' · 标准版'), note)
    w = a['w_rad_s']
    for axis, col in zip(ax, ['asym_mag_db', 'asym_phase_deg']):
        line, = axis.plot(w, a[col], color='#42627c' if sketch else BLUE, lw=1.9 if sketch else 2.1)
        if sketch:
            line.set_sketch_params(scale=.55, length=70, randomness=1.4)
        for wc in (2, 5):
            axis.axvline(wc, color=GRAY, ls=':', lw=.9)
    xmag = np.array([.2, 2., 5.])
    ymag = -20*np.log10(xmag)-20*np.log10(np.maximum(1, xmag/2))-20*np.log10(np.maximum(1, xmag/5))
    ax[0].scatter(xmag, ymag, color=BLUE, s=32, zorder=5)
    for x, y, offset in zip(xmag, ymag, [(10, 10), (-65, 15), (12, 7)]):
        ax[0].annotate(f'{x:g} rad/s\n{y:.2f} dB', (x, y), xytext=offset,
                       textcoords='offset points', fontsize=9,
                       arrowprops={'arrowstyle': '-', 'color': GRAY, 'lw': .6})
    for x, slope, dy in [(.46, '−20 dB/dec', 8), (3.1, '−40 dB/dec', -17), (20, '−60 dB/dec', 12)]:
        y = np.interp(np.log10(x), np.log10(w), a['asym_mag_db'])
        ax[0].annotate(slope, (x, y), xytext=(0, dy), textcoords='offset points',
                       ha='center', fontsize=10, color=BLUE)
    # All phase vertices are fixed markers, including overlapping transitions.
    xp = np.array([.1, .2, .5, 2., 5., 20., 50., 100.])
    yp = -90-45*np.clip(np.log10(xp/2)+1, 0, 2)-45*np.clip(np.log10(xp/5)+1, 0, 2)
    ax[1].scatter(xp, yp, color=BLUE, s=23, zorder=5)
    ax[1].annotate('相位变化区间：0.2—20 与 0.5—50 rad/s',
                   (.5, -104), xytext=(12, -2), textcoords='offset points', fontsize=9)
    for x in (.2, .5, 20, 50):
        ax[1].axvline(x, color=GRAY, ls=':', lw=.55, alpha=.7)
    ax[0].set_ylim(-110, 30)
    ax[0].yaxis.set_major_locator(MultipleLocator(20))
    ax[1].set_ylim(-285, -75)
    ax[1].yaxis.set_major_locator(MultipleLocator(45))
    save(fig, '01_asymptotes_' + ('sketch' if sketch else 'standard'))


def factors(a):
    fig, ax = pair('任务 A｜典型因子的幅值与相位贡献',
                   '归一化分解：G(s) = (1/s) × 1/(1+s/2) × 1/(1+s/5)；dB 与相位分别相加。')
    w = a['w_rad_s']
    for prefix, color, label in [('integrator', BLUE, r'$1/s$'),
                                  ('pole2', ORANGE, r'$1/(1+s/2)$'),
                                  ('pole5', GREEN, r'$1/(1+s/5)$')]:
        ax[0].plot(w, a[prefix+'_mag_db'], color=color, lw=1.6, label=label)
        ax[1].plot(w, a[prefix+'_phase_deg'], color=color, lw=1.6)
    ax[0].plot(w, a['exact_mag_db'], color='#252a30', lw=2, ls='--', label='因子之和')
    ax[1].plot(w, a['exact_phase_deg'], color='#252a30', lw=2, ls='--')
    ax[0].legend(ncol=4, loc='lower left', fontsize=9)
    save(fig, '02_factor_contributions')


def exact_bode(b, p):
    fig, ax = pair('任务 B｜精确 Bode 曲线与六个指定频率',
                   '圆点为实际计算读数；虚线给出同一模型的渐近近似。')
    for axis, field in zip(ax, ['mag_db', 'phase_deg']):
        axis.plot(b['w_rad_s'], b['exact_'+field], color=BLUE, lw=2, label='精确曲线')
        axis.plot(b['w_rad_s'], b['asym_'+field], color=ORANGE, ls='--', lw=1.5, label='渐近线')
        axis.scatter(p['w_rad_s'], p['exact_'+field], color=BLUE, edgecolor='white', s=42, zorder=5)
        for idx, (w, y) in enumerate(zip(p['w_rad_s'], p['exact_'+field])):
            offset = [(4, 12), (-42, 20), (-12, -32), (12, 10), (12, 15), (-65, -24)][idx]
            unit = 'dB' if field == 'mag_db' else '°'
            axis.annotate(f'{w:g}: {y:.2f} {unit}', (w, y), textcoords='offset points',
                          xytext=offset, fontsize=8.5,
                          arrowprops={'arrowstyle': '-', 'color': GRAY, 'lw': .6})
    ax[0].legend(loc='lower left', ncol=2)
    ax[0].set_ylim(-100, 30)
    ax[1].set_ylim(-290, -75)
    save(fig, '03_exact_bode_six_points')


def approximation_errors(b, p):
    fig, ax = pair('任务 B｜渐近线相对于精确曲线的误差',
                   '误差定义：渐近值 − 精确值；极值取实际频率网格上的最大绝对误差。')
    for axis, col, ylabel in zip(ax, ['mag_error_db', 'phase_error_deg'], ['幅值误差 (dB)', '相位误差 (°)']):
        axis.set_ylabel(ylabel)
        axis.plot(b['w_rad_s'], b[col], color=BLUE, lw=1.8)
        axis.scatter(p['w_rad_s'], p[col], color=ORANGE, s=26, zorder=5, label='六个指定频率')
        maximum = np.max(np.abs(b[col]))
        tied = np.flatnonzero(np.isclose(np.abs(b[col]), maximum, rtol=0, atol=1e-9))
        for k in tied:
            w, v = b['w_rad_s'][k], b[col][k]
            axis.scatter(w, v, marker='*', s=115, color=GREEN, zorder=6)
            offset = (-95, -25) if w > 20 else ((-85, 25) if w < 3 and w > 1 else (12, 25))
            axis.annotate(f'ω = {w:.4g} rad/s\n误差 = {v:+.4f}',
                          (w, v), xytext=offset, textcoords='offset points',
                          fontsize=9, arrowprops={'arrowstyle': '->', 'color': GREEN})
        axis.text(.98, .96, f'最大 |误差| = {maximum:.4f}'+(' dB' if col == 'mag_error_db' else '°'),
                  transform=axis.transAxes, ha='right', va='top', color=GREEN, fontsize=10)
        axis.axhline(0, color=GRAY, lw=.8)
        for wc in (2, 5):
            axis.axvline(wc, color=GRAY, ls=':', lw=.8)
        axis.margins(y=.35)
    ax[0].legend(loc='upper left')
    save(fig, '04_asymptote_errors')


MODELS = [('coarse', ORANGE, '斜率粗估', '--'),
          ('full16', GREEN, '16 点全精度拟合', '-.'),
          ('rounded16', '#7c55a4', '16 点舍入数据拟合', ':'),
          ('rounded8', '#b85a68', '8 点舍入数据拟合', '--')]


def identification(c, inputs):
    fig, ax = pair('任务 C｜离散频响辨识与模型验证',
                   '散点为 16 个全精度输入频率点；模型曲线由离散数据辨识后在稠密网格验证。')
    for axis, field in zip(ax, ['mag_db', 'phase_deg']):
        axis.plot(c['w_rad_s'], c['true_'+field], color=BLUE, lw=2.2, label='原模型')
        for key, color, label, ls in MODELS[:3]:
            axis.plot(c['w_rad_s'], c[key+'_'+field], color=color, ls=ls, lw=1.5, label=label)
        axis.scatter(inputs['w_rad_s'], inputs[field], facecolors='white', edgecolors=BLUE,
                     s=32, zorder=6, label='16 点输入' if field == 'mag_db' else None)
    ax[0].legend(loc='lower left', fontsize=8.5, ncol=2)
    save(fig, '05_identification_bode')

    fig, ax = pair('任务 C｜估计模型的幅相残差',
                   '残差定义：估计模型 − 原模型；全精度输入用于检验数值恢复，舍入输入用于检验有限精度影响。')
    for axis, field, ylabel in zip(ax, ['mag_db', 'phase_deg'], ['幅值残差 (dB)', '相位残差 (°)']):
        axis.set_ylabel(ylabel)
        for key, color, label, ls in MODELS:
            delta = c[key+'_'+field]-c['true_'+field]
            axis.plot(c['w_rad_s'], delta, color=color, ls=ls, lw=1.6,
                      label=f'{label}｜最大 |误差| {np.max(np.abs(delta)):.3g}')
        axis.axhline(0, color=GRAY, lw=.7)
        axis.legend(loc='best', fontsize=8)
    save(fig, '06_identification_residuals')


def sampling(c, full, rounded, sparse):
    fig, axs = plt.subplots(2, 2, figsize=(12.2, 8.4))
    fig.subplots_adjust(left=.085, right=.965, top=.865, bottom=.12, hspace=.34, wspace=.26)
    fig.suptitle('扩展验证｜离散采样数量与数据精度的影响', fontsize=15, y=.965)
    fig.text(.5, .915, '同一低阶结构，对比 16 点全精度、16 点舍入及 8 点舍入输入', ha='center', color=GRAY)
    for row, field, unit in [(0, 'mag_db', 'dB'), (1, 'phase_deg', '°')]:
        ax = axs[row, 0]
        decorate(ax, f'输入舍入误差 ({unit})')
        delta = rounded[field]-full[field]
        ax.plot(full['w_rad_s'], delta, 'o-', color='#7c55a4', ms=4, lw=1, label='16 点舍入 − 全精度')
        # Eight-point errors evaluated against full-precision points at identical frequencies.
        true_sparse = np.interp(np.log10(sparse['w_rad_s']), np.log10(full['w_rad_s']), full[field])
        # The numerical producer uses an actual subset; enforce exact membership.
        if not all(np.any(np.isclose(w, full['w_rad_s'], rtol=1e-12)) for w in sparse['w_rad_s']):
            raise ValueError('8-point input must be a subset of full16 for paired rounding errors')
        ax.scatter(sparse['w_rad_s'], sparse[field]-true_sparse, marker='x', color='#b85a68', s=40, label='8 点舍入 − 全精度')
        ax.axhline(0, color=GRAY, lw=.7)
        ax.legend(fontsize=8, loc='best')
        ax.set_title('输入数据的有限精度', fontsize=11)
        ax = axs[row, 1]
        decorate(ax, f'拟合残差 ({unit})')
        for key, color, label, ls in MODELS[1:]:
            delta = c[key+'_'+field]-c['true_'+field]
            ax.plot(c['w_rad_s'], delta, color=color, ls=ls, label=label, lw=1.5)
        ax.axhline(0, color=GRAY, lw=.7)
        ax.legend(fontsize=8, loc='best')
        ax.set_title('估计模型在连续频段上的误差', fontsize=11)
    fig.text(.085, .035, '扩展证据：舍入与减采样分别改变数据约束；残差由实际拟合模型在 0.1—100 rad/s 上计算。', fontsize=9, color=GRAY)
    save(fig, '07_sampling_precision')


def main():
    a, b, p = [load(n) for n in ['a_factors', 'b_dense', 'b_points']]
    c, full, rounded, sparse = [load(n) for n in ['c_validation', 'c_input_full16', 'c_input_rounded16', 'c_input_rounded8']]
    np.testing.assert_allclose(a['integrator_mag_db']+a['pole2_mag_db']+a['pole5_mag_db'], a['exact_mag_db'], atol=1e-9)
    np.testing.assert_allclose(a['integrator_phase_deg']+a['pole2_phase_deg']+a['pole5_phase_deg'], a['exact_phase_deg'], atol=1e-9)
    np.testing.assert_allclose(p['w_rad_s'], [.2, 1, 2, 5, 10, 50], atol=1e-12)
    asymptotes(a)
    asymptotes(a, sketch=True)
    factors(a)
    exact_bode(b, p)
    approximation_errors(b, p)
    identification(c, full)
    sampling(c, full, rounded, sparse)


if __name__ == '__main__':
    main()
