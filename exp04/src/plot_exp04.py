"""Render Octave-produced CSV evidence. Requires NumPy and Matplotlib.
Run this file with Python; all paths are derived from its location.
"""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.local' / 'exp04'
os.environ['MPLCONFIGDIR'] = str(PRIVATE / 'mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = PRIVATE / 'results'
plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                     'axes.spines.right': False, 'svg.fonttype': 'none'})
colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7']
labels = ['K = 1 (stable)', 'K = 5.9 (stable)', 'K = 6 (marginal)', 'K = 8 (unstable)']
def read(name):
    return np.loadtxt(OUT / name, delimiter=',', skiprows=1)
def save(fig, name):
    fig.savefig(OUT / (name+'.png'), dpi=180)
    fig.savefig(OUT / (name+'.svg'))
    plt.close(fig)
def decorate(ax, xlabel, ylabel):
    ax.set(xlabel=xlabel, ylabel=ylabel)
    ax.grid(alpha=.22)

p = read('poles.csv')
fig, axs = plt.subplots(1, 2, figsize=(12, 5), layout='constrained')
for ax in axs:
    ax.axvline(0, color='black', lw=1)
    ax.axhline(0, color='gray', lw=.7)
    for j, k in enumerate([1, 5.9, 6, 8]):
        a = p[p[:,0] == k]
        ax.scatter(a[:,1], a[:,2], marker=['o','s','x','^'][j], s=65,
                   color=colors[j], label=labels[j], zorder=4)
    decorate(ax, 'Real part (1/s)', 'Imaginary part (rad/s)')
axs[0].set_title('All closed-loop poles')
axs[0].legend(loc='upper left', fontsize=9)
axs[1].set(xlim=(-.025,.10), ylim=(1.38,1.61), title='Near the imaginary axis: upper pole')
axs[1].annotate('K = 6: +j sqrt(2)', (0,np.sqrt(2)), xytext=(.025,1.40),
                arrowprops={'arrowstyle':'->'}, fontsize=10)
fig.suptitle('Experiment 04 | Pole evidence and stability boundary')
save(fig,'01_poles')

a=read('step_overlay.csv'); long=read('step_long.csv')
fig, axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for j in range(4): axs[0,0].plot(a[:,0],a[:,j+1],color=colors[j],label=labels[j])
axs[0,0].set(xlim=(0,20),ylim=(-4,6),title='All four gains: common early-time window')
axs[0,0].legend(fontsize=9,loc='upper left')
for j in range(3): axs[0,1].plot(long[:,0],long[:,j+1],color=colors[j],label=labels[j],lw=.9)
axs[0,1].set(xlim=(0,60),title='Stable and marginal responses, 0-60 s')
for j in [1,2]: axs[1,0].plot(long[:,0],long[:,j+1],color=colors[j],label=labels[j],lw=.9)
axs[1,0].set(xlim=(1100,1200),title='Late-time distinction: decay vs persistent oscillation')
axs[1,0].legend(fontsize=9)
axs[1,1].plot(a[:,0],a[:,4],color=colors[3])
axs[1,1].set(title='K = 8: growing oscillation',xlim=(0,60))
for ax in axs.flat:
    ax.axhline(1,color='gray',ls='--',lw=.8)
    decorate(ax,'Time (s)','Output (unit step)')
fig.suptitle('Experiment 04 | Zero-state unit-step responses')
save(fig,'02_step_responses')

scan=read('scan_1.csv'); tracks=read('pole_trajectories.csv')
fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
axs[0].plot(scan[:,0],scan[:,1],color=colors[0],label='roots, step 0.05')
axs[0].axhline(0,color='black',lw=.8); axs[0].axvline(6,color=colors[2],ls='--',label='Routh boundary K = 6')
axs[0].legend(fontsize=9)
decorate(axs[0],'Gain K','Maximum pole real part (1/s)')
axs[0].set_title('Boundary crossing over K = 4...8')
sc=axs[1].scatter(tracks[:,1],tracks[:,2],c=tracks[:,0],s=15,cmap='viridis')
axs[1].axvline(0,color='black',lw=.8)
decorate(axs[1],'Real part (1/s)','Imaginary part (rad/s)')
axs[1].set_title('Parameter-pole trajectories')
fig.colorbar(sc,ax=axs[1],label='Gain K')
save(fig,'03_boundary_and_trajectories')

summary=read('scan_summary.csv'); z=read('K0_internal.csv')
fig,axs=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
axs[0].loglog(summary[:,0],summary[:,7],'o-',color=colors[0])
for row in summary:
    axs[0].annotate(f'[{row[3]:g}, {row[4]:g}]',(row[0],row[7]),xytext=(8,10 if row[0]<.001 else -16),textcoords='offset points',fontsize=9)
decorate(axs[0],'Scan increment in K','Half-width of sign bracket')
axs[0].set(xlim=(.0003,.15),title='Tenfold refinement reduces bracket uncertainty')
axs[1].plot(z[:,0],z[:,1],label='Initial y = 1, zero input',color=colors[1])
axs[1].plot(z[:,0],z[:,2],label='Unit step, zero initial state',color=colors[0])
axs[1].set(ylim=(-.2,1.3),title='K = 0: persistent internal mode and zero transfer')
axs[1].legend(fontsize=9)
decorate(axs[1],'Time (s)','Output')
save(fig,'04_refinement_and_internal_mode')
print('Saved four PNG/SVG figure pairs from Octave CSV data.')
