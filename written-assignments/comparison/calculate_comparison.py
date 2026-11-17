"""Independent algebra, numerical data and figures for two control problems."""
from pathlib import Path
import os
ROOT = Path(__file__).resolve().parent
os.environ['MPLCONFIGDIR'] = str(ROOT / '.local' / 'cache' / 'matplotlib')
import json
import numpy as np
import sympy as sp
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

x, v = sp.symbols('x v', real=True)
I = sp.I
line = -v + I*sp.Rational(4,3)*v

def cpair(z):
    return [float(complex(z).real), float(complex(z).imag)]

def coeff(expr):
    return np.array([float(a) for a in sp.Poly(expr, x).all_coeffs()])

def roots(expr):
    return np.roots(coeff(expr))

def state(expr,scale=1.):
    r = roots(expr)*scale
    zeta = -r.real / np.abs(r)
    return {'roots': [cpair(z) for z in sorted(r,key=lambda z:(z.real,z.imag))],
            'minimum_zeta':float(zeta.min()), 'stable':bool(np.all(r.real<0))}

def real_roots_exact(poly, symbol):
    # SymPy's exact real-root isolation proves all real candidates are counted.
    poly = sp.Poly(poly, symbol)
    return [float(sp.N(r,35)) for r in poly.real_roots()]

def damping_intersections(D,N,scale):
    ds = sp.expand(D.subs(x,line)); ns = sp.expand(N.subs(x,line))
    cross = sp.factor(sp.im(sp.expand(ds*sp.conjugate(ns))))
    raw = []
    for t in real_roots_exact(cross,v):
        if t <= 0: continue
        z = complex(line.subs(v,t))
        gain = -complex(D.subs(x,z))/complex(N.subs(x,z))
        if gain.real>0:
            raw.append({'sigma_normalized': t, 'sigma':t*scale,
                        'gain':gain.real,'imaginary_gain_residual':gain.imag,
                        'boundary_pole':[z.real*scale,z.imag*scale]})
    raw.sort(key=lambda a:a['gain'])
    return {'exact_cross_polynomial':str(cross),'positive_gain_crossings':raw}

def breakpoints(D,N,scale):
    poly = sp.expand(sp.diff(D,x)*N-D*sp.diff(N,x))
    result=[]
    for exact_t in sp.Poly(poly,x).real_roots():
        t=float(sp.N(exact_t,45))
        nv=float(sp.N(N.subs(x,exact_t),45))
        if abs(nv)<1e-25: continue
        gain=float(sp.N(-D.subs(x,exact_t)/N.subs(x,exact_t),45))
        if gain>0:
            result.append({'root':t*scale,'normalized_root':t,'gain':gain})
    return result

def imaginary_crossings(D,N,scale):
    z=I*v
    ds=sp.expand(D.subs(x,z)); ns=sp.expand(N.subs(x,z))
    cross=sp.factor(sp.im(sp.expand(ds*sp.conjugate(ns))))
    result=[]
    for w in real_roots_exact(cross,v):
        if w<=0:continue
        gain=-complex(D.subs(x,I*w))/complex(N.subs(x,I*w))
        result.append({'frequency':w*scale,'normalized_frequency':w,
                       'gain':gain.real,'imaginary_gain_residual':gain.imag})
    return {'exact_cross_polynomial':str(cross),'crossings':result}

def interval_tests(D,N,crossings,scale=1.):
    b=[0]+sorted(c['gain'] for c in crossings)+[float('inf')]
    out=[]
    for lo,hi in zip(b[:-1],b[1:]):
        g=(hi/2 if lo==0 else (lo*2 if np.isinf(hi) else np.sqrt(lo*hi)))
        item={'interval':[lo, None if np.isinf(hi) else hi], 'test_gain':g}
        item.update(state(D+g*N,scale));out.append(item)
    return out

def root_direction(D,N,pole):
    q=sp.diff(D,x).subs(x,pole)
    dz=-complex(N.subs(x,pole))/complex(q)
    return float(np.degrees(np.angle(dz)))

models={
 'seconds':{'e':sp.Rational(8,10**7),'u':sp.Rational(2,10**6),
            'b':sp.Rational(4,10**10),'scale':1000.,'Ka':1000.},
 'millisecond_parameter_assumption':{'e':sp.Rational(4,5),'u':sp.Integer(2),
            'b':sp.Rational(2,5),'scale':1.,'Ka':.02}}

report={'calculation':'Exact rational polynomial construction and real-root isolation; independent numerical roots.'}
for name,m in models.items():
    D=sp.expand(x*(x+1)**2*(x*x+2*x+m['e']))
    Na=sp.expand(10*m['u']*x*(x+1)+m['b'])
    ka=sp.Rational(str(m['Ka']))
    Dc=sp.expand(D+m['b']*ka)
    Nc=sp.expand(m['u']*ka*x*(x+1))
    m.update(D=D,Na=Na,Dc=Dc,Nc=Nc)
    item={'normalized_characteristic':str(D)+' + Ka*('+str(Na)+')',
          'root_unit':'s^-1' if name=='seconds' else 'ms^-1',
          'root_multiplier_to_declared_time_variable':m['scale'],
          'Ka_open_poles':[cpair(z*m['scale']) for z in roots(D)],
          'Ka_open_zeros':[cpair(z*m['scale']) for z in roots(Na)],
          'Ka_breakpoints':breakpoints(D,Na,m['scale']),
          'Ka_imaginary':imaginary_crossings(D,Na,m['scale']),
          'Ka_damping':damping_intersections(D,Na,m['scale']),
          'chosen_Ka':m['Ka'], 'chosen_K2':10,
          'chosen_state':state(D+ka*Na,m['scale']),
          'K2_open_poles':[cpair(z*m['scale']) for z in roots(Dc)],
          'K2_open_zeros':[cpair(z*m['scale']) for z in roots(Nc)],
          'K2_breakpoints':breakpoints(Dc,Nc,m['scale']),
          'K2_imaginary':imaginary_crossings(Dc,Nc,m['scale']),
          'K2_damping':damping_intersections(Dc,Nc,m['scale'])}
    for axis,dd,nn in [('Ka',D,Na),('K2',Dc,Nc)]:
        cs=item[axis+'_damping']['positive_gain_crossings']
        item[axis+'_damping_interval_tests']=interval_tests(dd,nn,cs,m['scale'])
    item['K2_initial_departure_angles']= [
        {'pole':cpair(z*m['scale']), 'angle_degrees':root_direction(Dc,Nc,z)}
        for z in roots(Dc) if z.imag>1e-8]
    report[name]=item

# Verify the time-unit transformation symbolically for arbitrary gains.
s, Ka, K2 = sp.symbols('s Ka K2', real=True)
Qs=s*(s/1000+1)**2*(s*s+2000*s+sp.Rational(4,5))+2*Ka*K2*s*(s/1000+1)+sp.Rational(2,5)*Ka
expected=x*(x+1)**2*(x*x+2*x+sp.Rational(8,10**7))+sp.Rational(2,10**6)*Ka*K2*x*(x+1)+sp.Rational(4,10**10)*Ka
report['unit_scaling_identity_exact']= bool(sp.expand(Qs.subs(s,1000*x)/10**9-expected)==0)
report['unit_scaling_polynomial']=str(sp.expand(expected))

# AP6.5: an exact example and independent step response from residues.
Da=2*x**3+64*x*x+42*x+180
ra=roots(Da)
ts=np.linspace(0,35,7001)
t,ys=signal.step(signal.TransferFunction([30],coeff(Da)),T=ts)
res,pole,_=signal.residue([30],np.r_[coeff(Da),0.])
yr=(np.exp(np.outer(t,pole))@res).real
err=float(np.max(np.abs(ys-yr)))
assert err<1e-10
np.savetxt(ROOT/'AP6_5_step.csv',np.c_[t,ys,yr],delimiter=',',header='time_s,step_state_space,step_residue',comments='')
report['AP6_5']={'gains_K1_to_K5':[1,1,20,300,46],
                'a':21,'b':42,'g':9000,'D':[2,64,42,180],
                'poles':[cpair(z) for z in ra],'final_value':1/6,
                'step_independent_max_difference':err,
                'step_peak_time_s':float(t[np.argmax(ys)]),
                'step_peak':float(np.max(ys))}

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.grid':True,
                     'grid.alpha':.25,'savefig.dpi':180})
fig,ax=plt.subplots(figsize=(8,4.3),constrained_layout=True)
ax.plot(t,ys,label='Unit-step response')
ax.axhline(1/6,color='0.35',ls='--',label='Final value 1/6')
ax.set(xlabel='Time (s)',ylabel='Output',title='AP6.5(a): K1=1, K2=1, K3=20, K4=300, K5=46')
ax.legend(); fig.savefig(ROOT/'AP6_5_step.png');fig.savefig(ROOT/'AP6_5_step.svg');plt.close(fig)

# Publish all sampled polynomial roots and damping minima, with no external input.
for name,m in models.items():
    for axis,dd,nn,limits in [('Ka',m['D'],m['Na'],(-9,9) if name=='seconds' else (-7,3)),
                              ('K2',m['Dc'],m['Nc'],(-8,7) if name=='seconds' else (-5,5))]:
        gains=np.r_[0.,np.logspace(*limits,2200)]
        rr=np.array([roots(dd+float(g)*nn)*m['scale'] for g in gains])
        zz=np.min(-rr.real/np.maximum(np.abs(rr),np.finfo(float).tiny),axis=1)
        tab=np.array([[g,j,z.real,z.imag,zz[i]] for i,(g,rs) in enumerate(zip(gains,rr)) for j,z in enumerate(rs)])
        np.savetxt(ROOT/(name+'_'+axis+'_locus.csv'),tab,delimiter=',',
                   header='gain,root_index,root_real,root_imag,minimum_zeta',comments='')
        m[axis+'_gains']=gains;m[axis+'_roots']=rr;m[axis+'_zeta']=zz

for axis in ['Ka','K2']:
    fig,axs=plt.subplots(2,2,figsize=(11,8),constrained_layout=True)
    spec=[('seconds',axs[0,0],(-2400,300),(-1700,1700),'Seconds: full view'),
          ('seconds',axs[0,1],(-1.1,.3),(-.8,.8),'Seconds: slow modes'),
          ('millisecond_parameter_assumption',axs[1,0],(-3.5,1),(-2.5,2.5),'Separate millisecond parameter assumption')]
    for name,ax,xlim,ylim,title in spec:
        m=models[name]; rr=m[axis+'_roots']
        ax.plot(rr.real.ravel(),rr.imag.ravel(),'.',color='#2079ac',ms=1.1,alpha=.8)
        pp=roots(m['D'] if axis=='Ka' else m['Dc'])*m['scale']
        zz=roots(m['Na'] if axis=='Ka' else m['Nc'])*m['scale']
        ax.plot(pp.real,pp.imag,'x',color='#be391d',ms=7,mew=1.4)
        ax.plot(zz.real,zz.imag,'o',color='#267838',mfc='none',ms=6)
        u=np.linspace(xlim[0],0,100)
        ax.plot(u,-4*u/3,'--',color='0.55',lw=.7)
        ax.plot(u,4*u/3,'--',color='0.55',lw=.7)
        ax.axvline(0,color='0.4',lw=.6);ax.axhline(0,color='0.4',lw=.6)
        ax.set(xlim=xlim,ylim=ylim,title=title,xlabel='Real part ('+('s^-1' if name=='seconds' else 'ms^-1')+')',ylabel='Imaginary part')
    ax=axs[1,1]
    for name,label,color in [('seconds','Seconds','#2079ac'),('millisecond_parameter_assumption','Separate ms assumption','#bc5826')]:
        m=models[name];ax.semilogx(m[axis+'_gains'][1:],m[axis+'_zeta'][1:],label=label,color=color)
    ax.axhline(.6,color='0.25',ls='--',lw=1,label='Required damping 0.6')
    ax.set(xlabel=axis,ylabel='Minimum pole damping ratio',ylim=(-.6,1.05),title='All five poles evaluated')
    ax.legend(loc='lower left')
    fig.savefig(ROOT/('P7_11_'+axis+'_comparison.png'));fig.savefig(ROOT/('P7_11_'+axis+'_comparison.svg'));plt.close(fig)

(ROOT/'numerical_results.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({name:{key:report[name][key] for key in ['Ka_breakpoints','Ka_damping','Ka_imaginary','K2_breakpoints','K2_imaginary','K2_damping']} for name in models},indent=2))
if (ROOT/'plot_existing.py').is_file():
    import runpy
    runpy.run_path(str(ROOT/'plot_existing.py'),run_name='__main__')
