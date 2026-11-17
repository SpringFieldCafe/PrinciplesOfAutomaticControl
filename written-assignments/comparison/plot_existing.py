"""Redraw existing numerical data; arrows follow increasing gain."""
from pathlib import Path
import os
ROOT=Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR']=str(ROOT/'.local'/'cache'/'matplotlib')
import json
import numpy as np
from scipy.optimize import linear_sum_assignment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

data=json.loads((ROOT/'numerical_results.json').read_text(encoding='utf-8'))
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,
                     'grid.alpha':.25,'savefig.dpi':180})
names=['seconds','millisecond_parameter_assumption']
for axis in ['Ka','K2']:
    fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    paths={}
    for name in names:
        tab=np.loadtxt(ROOT/(name+'_'+axis+'_locus.csv'),delimiter=',',skiprows=1)
        gain=tab[::5,0]
        rr=(tab[:,2]+1j*tab[:,3]).reshape(-1,5)
        branches=np.empty_like(rr);branches[0]=rr[0]
        for k in range(1,len(gain)):
            rows,cols=linear_sum_assignment(abs(branches[k-1,:,None]-rr[k,None,:]))
            branches[k]=rr[k,cols]
        paths[name]=(gain,branches,tab[::5,4])
    ka=axis=='Ka'
    spec=[('seconds',axs[0,0],(-2400,300),(-1700,1700),'Seconds: full view',
           [1e3,1e5,1e6] if ka else [100,1000,10000]),
          ('seconds',axs[0,1],(-.06,.005) if ka else (-1.1,.3),(-.025,.025) if ka else (-.8,.8),
           'Seconds: slow modes',[.01,2,100] if ka else [.1,1,5]),
          ('millisecond_parameter_assumption',axs[1,0],(-3.5,1),(-2.5,2.5),
           'Separate millisecond parameter assumption',[.01,.1,1] if ka else [10,100,300])]
    for name,ax,xlim,ylim,title,targets in spec:
        gain,rr,_=paths[name]
        ax.plot(rr.real,rr.imag,'.',color='#2079ac',ms=1.2)
        for target in targets:
            ia=int(np.argmin(abs(gain-target)))
            ib=int(np.argmin(abs(gain-target*1.55)))
            for j in range(5):
                a,b=rr[ia,j],rr[ib,j]
                visible=all(xlim[0]<z.real<xlim[1] and ylim[0]<z.imag<ylim[1] for z in [a,b])
                length=np.hypot((a.real-b.real)/(xlim[1]-xlim[0]),(a.imag-b.imag)/(ylim[1]-ylim[0]))
                if visible and length>.006:
                    ax.annotate('',xy=(b.real,b.imag),xytext=(a.real,a.imag),
                                arrowprops={'arrowstyle':'->','color':'#164c75','lw':1.2,'mutation_scale':10})
        item=data[name]
        pp=np.array(item['Ka_open_poles' if ka else 'K2_open_poles'])
        zz=np.array(item['Ka_open_zeros' if ka else 'K2_open_zeros'])
        ax.plot(pp[:,0],pp[:,1],'x',color='#be391d',ms=7,mew=1.4)
        ax.plot(zz[:,0],zz[:,1],'o',color='#267838',mfc='none',ms=6)
        u=np.linspace(xlim[0],0,100)
        ax.plot(u,-4*u/3,'--',color='0.55',lw=.7)
        ax.plot(u,4*u/3,'--',color='0.55',lw=.7)
        ax.axvline(0,color='0.4',lw=.6);ax.axhline(0,color='0.4',lw=.6)
        unit=r'$\mathrm{s}^{-1}$' if name=='seconds' else r'$\mathrm{ms}^{-1}$'
        ax.set(xlim=xlim,ylim=ylim,title=title,xlabel='Real part ('+unit+')',ylabel='Imaginary part ('+unit+')')
    ax=axs[1,1]
    for name,label,color in [('seconds','Seconds','#2079ac'),('millisecond_parameter_assumption','Separate ms assumption','#bc5826')]:
        gain,_,zeta=paths[name]
        ax.semilogx(gain[1:],zeta[1:],label=label,color=color)
    ax.axhline(.6,color='0.25',ls='--',lw=1,label='Required damping 0.6')
    ax.set(xlabel=r'$K_a$' if ka else r'$K_2$',ylabel='Minimum pole damping ratio',ylim=(-.6,1.05),title='All five poles evaluated')
    ax.legend(loc='lower left')
    fig.savefig(ROOT/('P7_11_'+axis+'_comparison.png'))
    fig.savefig(ROOT/('P7_11_'+axis+'_comparison.svg'))
    plt.close(fig)
print('Redrew two comparison figures from existing CSV data.')
