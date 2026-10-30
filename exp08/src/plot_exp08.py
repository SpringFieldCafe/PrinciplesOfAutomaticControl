"""Plot Octave CSV evidence and an editable SVG signal diagram."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'.local'/'exp08'/'results'
os.environ['MPLCONFIGDIR']=str(OUT.parent/'mplconfig')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
colors=['#777777','#0072B2','#E69F00','#009E73']
labels=['Baseline C=10','Lead A (selected)','Lead B (Kv=0.8)','Lead-lag (Kv=1.2)']
def load(n):return np.loadtxt(OUT/(n+'.csv'),delimiter=',',skiprows=1)
def save(fig,n):
    for ext in ['png','svg']:fig.savefig(OUT/(n+'.'+ext),dpi=180)
    plt.close(fig)
def deco(ax,x,y):ax.set(xlabel=x,ylabel=y);ax.grid(alpha=.22)
b=load('bode');m=load('metrics');y=load('step');st=load('sensitivity');r=load('ramp_error')
fig,axs=plt.subplots(2,1,figsize=(10,7),layout='constrained',sharex=True)
for j in range(4):
    axs[0].semilogx(b[:,0],b[:,1+2*j],label=labels[j],color=colors[j])
    axs[1].semilogx(b[:,0],b[:,2+2*j],color=colors[j])
    axs[0].scatter(m[j,6],0,color=colors[j],s=28,zorder=5)
axs[0].axhline(0,color='black',ls='--',lw=.8);axs[1].axhline(-180,color='black',ls='--',lw=.8)
axs[0].legend(fontsize=9);axs[0].set(ylim=(-80,60),title='Experiment 08 | Open-loop Bode comparison')
axs[1].set(xlim=(.001,100),ylim=(-280,-60))
deco(axs[0],'','Magnitude (dB)');deco(axs[1],'Angular frequency (rad/s)','Phase (deg)')
save(fig,'01_bode')
fig,axs=plt.subplots(1,2,figsize=(12,5),layout='constrained')
for j in range(4):
    axs[0].plot(y[:,0],y[:,j+1],color=colors[j],label=labels[j]);axs[1].plot(y[:,0],y[:,j+1],color=colors[j])
for ax in axs:
    ax.axhspan(.98,1.02,color='gray',alpha=.12);ax.axvline(4,color='black',ls='--',lw=.8)
    deco(ax,'Time (s)','Output (unit step)')
axs[0].set(xlim=(0,12),ylim=(0,1.18),title='Rise, overshoot and settling');axs[0].legend(fontsize=9)
axs[1].set(xlim=(0,30),ylim=(.97,1.035),title='2% band and lag tail')
save(fig,'02_step')
fig,axs=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
for j in range(4):
    axs[0].plot(r[:,0],r[:,j+1],color=colors[j],label=labels[j]);axs[0].axhline(1/m[j,1],color=colors[j],ls=':',lw=.7)
    axs[1].semilogx(st[:,0],20*np.log10(st[:,1+3*j]),color=colors[j])
    axs[2].semilogx(st[:,0],20*np.log10(st[:,3+3*j]),color=colors[j])
axs[0].set(xlim=(0,120),ylim=(0,1.4),title='Unit-ramp error: dotted = 1/Kv');axs[0].legend(fontsize=8)
axs[1].set(xlim=(.001,100),ylim=(-70,10),title='Sensitivity S: disturbance rejection')
axs[2].set(xlim=(.001,100),ylim=(-40,40),title='Control sensitivity CS: noise effort')
deco(axs[0],'Time (s)','Tracking error');deco(axs[1],'Angular frequency (rad/s)','Magnitude (dB)');deco(axs[2],'Angular frequency (rad/s)','Magnitude (dB)')
save(fig,'03_precision_and_sensitivity')
rob=load('robustness');fig,axs=plt.subplots(1,2,figsize=(11,5),layout='constrained')
for ax,id in zip(axs,[1,3]):
    a=rob[rob[:,0]==id];passed=a[:,10].astype(bool)
    ax.scatter(a[~passed,8],a[~passed,6],c='#D55E00',marker='x',s=65,label='One or more targets missed')
    ax.scatter(a[passed,8],a[passed,6],c='#009E73',marker='o',s=40,label='All four targets met')
    ax.axvline(4,color='black',ls='--');ax.axhline(1,color='black',ls='--')
    ax.set_title(f'{labels[id]}: {sum(passed)}/27 meet targets')
    deco(ax,'2% settling time (s)','Velocity constant Kv')
axs[0].legend(fontsize=8);fig.suptitle('Independent plant gain and pole changes: -20%, nominal, +20%')
save(fig,'04_robustness')

# Editable diagram; simulation counterpart is the C/G/feedback interconnection in run_exp08.m.
fig,ax=plt.subplots(figsize=(12,4),layout='constrained');ax.set(xlim=(0,12),ylim=(0,4));ax.axis('off')
def arrow(a,b):ax.annotate('',xy=b,xytext=a,arrowprops={'arrowstyle':'->','lw':1.5})
for x,text in [(3,'Controller C(s)\nselected: lead A'),(7,'Plant G(s)\n1/[s(s+2)(s+5)]')]:
    ax.add_patch(Rectangle((x,2),2.5,1,fill=False,lw=1.5));ax.text(x+1.25,2.5,text,ha='center',va='center')
ax.add_patch(Circle((1.7,2.5),.23,fill=False));ax.text(1.7,2.5,'+',ha='center',va='center');ax.text(1.35,2.1,'-',fontsize=15)
arrow((.3,2.5),(1.45,2.5));arrow((1.95,2.5),(3,2.5));arrow((5.5,2.5),(7,2.5));arrow((9.5,2.5),(11.5,2.5))
ax.plot([10.5,10.5,1.7],[2.5,.8,.8],color='black',lw=1.5);arrow((1.7,.8),(1.7,2.25))
for x,text in [(.6,'r'),(2.35,'e'),(6.1,'u'),(10.8,'y')]:ax.text(x,2.8,text)
ax.text(6,.45,'Unity negative feedback | T = CG/(1+CG), S = 1/(1+CG)',ha='center')
ax.set_title('Experiment 08 | Editable control signal diagram')
save(fig,'05_signal_diagram')
fig,ax=plt.subplots(figsize=(9,4.5),layout='constrained')
for j in range(4):ax.semilogx(st[:,0],20*np.log10(st[:,2+3*j]),color=colors[j],label=labels[j])
ax.set(xlim=(.01,100),ylim=(-80,10),title='Complementary sensitivity T: tracking bandwidth and output noise')
deco(ax,'Angular frequency (rad/s)','Magnitude (dB)');ax.legend(fontsize=9);save(fig,'06_complementary_sensitivity')
print('Saved six annotated PNG/SVG figure pairs.')
