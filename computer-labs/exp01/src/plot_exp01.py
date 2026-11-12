"""Plot actual experiment 01 data. NumPy + Matplotlib, no simulation here."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parents[2]
PRIVATE = ROOT / '.local' / 'exp01'
os.environ['MPLCONFIGDIR'] = str(PRIVATE / 'mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = PRIVATE / 'results'
plt.rcParams.update({'font.size': 11, 'svg.fonttype': 'none'})
COLORS = ['#0072B2', '#D55E00', '#009E73']
def read(name):
    return np.loadtxt(OUT / name, delimiter=',', skiprows=1)
def save(fig, stem):
    fig.savefig(OUT / (stem + '.png'), dpi=180)
    fig.savefig(OUT / (stem + '.svg'))
    plt.close(fig)
def axes(ax, title, ylabel='Displacement (m)'):
    ax.set(title=title, xlabel='Time (s)', ylabel=ylabel)
    ax.grid(alpha=.25)

a = read('open_loop_steps.csv')
fig, ax = plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
for j, name in enumerate(['Transfer function','Zero-pole-gain','Physical state space']):
    ax[0].plot(a[:,0],a[:,j+1],label=name,lw=2,ls=['-','--',':'][j],color=COLORS[j])
ax[0].legend(); axes(ax[0],'A | Equivalent models, 1 N force step')
for j in range(2,a.shape[1]): ax[1].plot(a[:,0],a[:,j]-a[:,1],label=f'Model {j} - tf')
axes(ax[1],'Conversion residuals', 'Response difference (m)'); ax[1].legend(fontsize=8)
save(fig,'01_model_equivalence')

a=read('closed_loop_steps.csv'); m=read('closed_loop_metrics.csv')
fig, ax=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
for j,K in enumerate([1,4,10]):
    ax[0].plot(a[:,0],a[:,j+1],label=f'K={K}',color=COLORS[j])
    ax[0].axhline(m[j,5],color=COLORS[j],ls=':',lw=.8)
    ax[1].plot(a[:,0],a[:,j+1]/m[j,5],label=f'K={K}',color=COLORS[j])
for aa in ax: aa.set_xlim(0,10); aa.legend()
axes(ax[0],'B | Unity negative feedback, 1 m reference step')
axes(ax[1],'Normalized response exposes overshoot','Output / steady output')
ax[1].axhline(1,color='gray',ls='--',lw=.8)
save(fig,'02_gain_comparison')

fig,ax=plt.subplots(figsize=(6,5),layout='constrained')
for j,K in enumerate([1,4,10]):
    ax.scatter(m[j,[1,3]],m[j,[2,4]],label=f'K={K}',color=COLORS[j],marker='x',s=90)
ax.axvline(0,color='black',lw=.8); ax.axhline(0,color='gray',lw=.8)
ax.set(xlabel='Real part (1/s)',ylabel='Imaginary part (rad/s)',title='B | Closed-loop poles: -1 +/- j sqrt(4+K)',xlim=(-1.5,.1))
ax.legend(); ax.grid(alpha=.25); save(fig,'03_closed_loop_poles')

# Xcos comparison is generated only after its actual CSV exports exist.
paths=list(OUT.glob('xcos_K*_h*.csv'))
if paths:
    fig,axs=plt.subplots(3,2,figsize=(12,10),layout='constrained')
    metrics=[]
    for j,K in enumerate([1,4,10]):
        for path in sorted(OUT.glob(f'xcos_K{K}_h*.csv')):
            x=np.loadtxt(path,delimiter=',',skiprows=1)
            h=float(path.stem.split('_h')[1]); t=x[:,0]; y=x[:,1]
            wd=np.sqrt(4+K); exact=K/(5+K)*(1-np.exp(-t)*(np.cos(wd*t)+np.sin(wd*t)/wd))
            octave=np.interp(t,a[:,0],a[:,j+1])
            assert np.max(np.abs(y-exact)) < 1e-5
            assert np.max(np.abs(y-octave)) < 1e-5
            assert np.max(np.abs(np.diff(t)-h)) < 1e-10
            metrics.append([K,h,len(t),t[-1],np.max(np.abs(y-exact)),np.max(np.abs(y-octave))])
            axs[j,0].plot(t,y,label=f'Xcos RK4, h={h:g} s')
            axs[j,1].plot(t,y-exact,label=f'h={h:g} s')
        axs[j,0].plot(a[:,0],a[:,j+1],'k--',lw=1,label='Octave step')
        axes(axs[j,0],f'C | K={K}: actual Xcos and Octave outputs')
        axes(axs[j,1],f'Xcos residual against analytic response, K={K}','Difference (m)')
        for ax in axs[j]: ax.legend(fontsize=9)
    np.savetxt(OUT/'xcos_errors.csv',metrics,delimiter=',',header='K,max_step_s,samples,last_recorded_time_s,max_abs_analytic_error_m,max_abs_octave_error_m',comments='')
    save(fig,'04_xcos_comparison')

# Clean signal diagram complements the native Xcos export.
from matplotlib.patches import Rectangle, Circle
fig,ax=plt.subplots(figsize=(12,4),layout='constrained')
ax.set(xlim=(0,12),ylim=(0,4)); ax.axis('off')
def arrow(start,end):
    ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','lw':1.6})
for x,width,label in [(3,1.6,'Gain K\n(N/m)'),(6,2.8,'Plant\n1 / (s² + 2s + 5)')]:
    ax.add_patch(Rectangle((x,2),width,1,fill=False,lw=1.6))
    ax.text(x+width/2,2.5,label,ha='center',va='center')
ax.add_patch(Circle((1.7,2.5),.22,fill=False,lw=1.6))
ax.text(1.7,2.5,'Σ',ha='center',va='center'); ax.text(1.3,2.75,'+')
ax.text(1.4,2.05,'−'); ax.text(.3,2.8,'r = 1 m step')
arrow((.3,2.5),(1.48,2.5)); arrow((1.92,2.5),(3,2.5))
arrow((4.6,2.5),(6,2.5)); arrow((8.8,2.5),(11.5,2.5))
ax.plot([10,10,1.7],[2.5,.8,.8],color='black',lw=1.6); arrow((1.7,.8),(1.7,2.28))
for x,label in [(2.3,'e (m)'),(5.1,'u (N)'),(10.5,'x (m)')]: ax.text(x,2.8,label)
ax.text(6,.35,'Zero initial state | Unity negative feedback | T(s) = K / (s² + 2s + 5 + K)',ha='center')
ax.set_title('Experiment 01 | Closed-loop signal relations')
save(fig,'06_signal_relations')
print('Experiment 01 figures generated from actual CSV files.')
