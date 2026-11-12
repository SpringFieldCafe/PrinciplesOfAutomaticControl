function run_exp05(make_plots)
% GNU Octave 11/control 4.2.1. Zero initial conditions; unit step; time s.
% Add exp05/src to Octave's path and call run_exp05. For short tool calls,
% run_exp05(false), then plot_exp05(1:2) and plot_exp05(3:5).
% All generated files stay inside the project, independently of startup pwd.
  if nargin<1, make_plots=true; end
  pkg load control;
  root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out=fullfile(root,'.local','exp05','results'); if ~exist(out,'dir'), mkdir(out); end
  diary(fullfile(out,'run.log')); diary on;
  printf('Octave %s; control package; unit step; zero state; dt=0.002 s\n',version); pkg list control;
  L0=tf(1,[1 6 8 0]);
  % Numerical noninteractive replacement for rlocfind on zeta=0.5 ray.
  gain_at=@(sigma) -polyval([1 6 8 0],-sigma+1i*sqrt(3)*sigma);
  sigma=fzero(@(x) imag(gain_at(x)),[.3 1]);
  selected=-sigma+1i*sqrt(3)*sigma; Ksel=real(gain_at(sigma));
  assert(abs(Ksel-224/27)<1e-10);
  assert(abs(1+Ksel/polyval([1 6 8 0],selected))<1e-12);
  candidates=roots([3 12 8]); gains=-polyval([1 6 8 0],candidates);
  writecsv(fullfile(out,'breakaway_candidates.csv'),'s_real,K,on_positive_gain_locus',[candidates,gains,gains>0]);
  [rdata,kbreak,rlpol,gvec]=rlocus(L0);
  residual=0;
  for j=1:numel(gvec), residual=max(residual,max(abs(polyval([1 6 8 gvec(j)],rlpol(:,j))))); end
  assert(residual<1e-9);
  writecsv(fullfile(out,'root_locus.csv'),'K,p1_real,p1_imag,p2_real,p2_imag,p3_real,p3_imag',[gvec(:),real(rlpol(1,:))',imag(rlpol(1,:))',real(rlpol(2,:))',imag(rlpol(2,:))',real(rlpol(3,:))',imag(rlpol(3,:))']);
  K=[5 Ksel 20 40 48 50]; names={'K5','selected','K20','K40','K48','K50'};
  times=cell(1,6); responses=cell(1,6); models=cell(1,6); metrics=NaN(6,4); predicted=NaN(6,2);
  mf=fopen(fullfile(out,'metrics.csv'),'w');
  fprintf(mf,'case,K,stable,dcgain,dominant_real,dominant_imag,zeta,wn_rad_s,Mp_percent,tp_s,tr10_90_s,ts2pct_s,Mp_pred_percent,ts_pred_s,Mp_actual_minus_pred_pp,ts_actual_minus_pred_s,horizon_s,final_error,modal_waveform_error,dt_Mp_delta_pp,dt_tp_delta_s,dt_tr_delta_s,dt_ts_delta_s\n');
  pf=fopen(fullfile(out,'poles.csv'),'w'); fprintf(pf,'case,K,pole_real,pole_imag,char_eq_residual\n');
  maxdiff=zeros(1,4); modalmax=0;
  for j=1:6
    kval=K(j); den=[1 6 8 kval]; p=roots(den); stable=all(real(p)<-1e-9);
    for h=1:3, fprintf(pf,'%s,%.15g,%.15g,%.15g,%.15g\n',names{j},kval,real(p(h)),imag(p(h)),abs(polyval(den,p(h)))); end
    sys=feedback(kval*L0,1); models{j}=sys;
    q=p(imag(p)>1e-8); assert(numel(q)==1); q=q(1); zeta=-real(q)/abs(q);
    if stable, horizon=max(40,ceil(24/(-real(q)))); else, horizon=80; end
    t=(0:.002:horizon)'; y=step(sys,t); times{j}=t; responses{j}=y;
    [rr,pp,kk]=residue(kval,conv(den,[1 0])); ya=zeros(size(t));
    for h=1:numel(pp), ya=ya+rr(h)*exp(pp(h)*t); end
    wave=max(abs(y-real(ya))); modalmax=max(modalmax,wave); assert(wave<1e-7);
    delta=NaN(1,4); finalerr=NaN;
    if stable
      yf=dcgain(sys); metrics(j,:)=measure(t,y,yf);
      tc=t(1:2:end); yc=step(sys,tc); delta=abs(measure(tc,yc,yf)-metrics(j,:));
      assert(delta(1)<.001 && max(delta(2:4))<1e-4); maxdiff=max(maxdiff,delta);
      finalerr=abs(y(end)-yf); assert(finalerr<1e-8);
      predicted(j,:)=[100*exp(-pi*zeta/sqrt(1-zeta^2)),4/(-real(q))];
    end
    m=metrics(j,:); pred=predicted(j,:);
    fprintf(mf,'%s',names{j}); fprintf(mf,',%.15g',[kval,stable,dcgain(sys),real(q),imag(q),zeta,abs(q),m,pred,m(1)-pred(1),m(4)-pred(2),horizon,finalerr,wave,delta]); fprintf(mf,'\n');
    writecsv(fullfile(out,[names{j} '_response.csv']),'time_s,output',[t,y]);
    printf('%s K=%.9g stable=%d: Mp=%.6f%% tp=%.6f tr=%.6f ts=%.6f s; poles real max=%.9g\n',names{j},kval,stable,m,max(real(p)));
  end
  fclose(mf); fclose(pf);
  assert(max(real(roots([1 6 8 47.99])))<0 && max(real(roots([1 6 8 48.01])))>0);
  boundary=roots([1 6 8 48]); assert(min(abs(boundary-1i*sqrt(8)))<1e-12);
  writecsv(fullfile(out,'selection.csv'),'K,selected_real,selected_imag,zeta,angle_condition_residual,characteristic_residual',[Ksel,real(selected),imag(selected),-real(selected)/abs(selected),abs(imag(gain_at(sigma))),abs(polyval([1 6 8 Ksel],selected))]);
  save(fullfile(out,'exp05_data.mat'),'times','responses','metrics','predicted','K','names','Ksel','selected','rlpol','gvec','candidates','gains');
  vf=fopen(fullfile(out,'validation.txt'),'w');
  fprintf(vf,'PASS: numerical zeta=.5 selection agrees with K=224/27; characteristic and angle conditions <1e-12.\nPASS: rlocus characteristic residual %.9g.\nPASS: independent partial-fraction step responses, max error %.9g.\nPASS: dt .004 -> .002 s max metric changes: %.9g pp, %.9g s, %.9g s, %.9g s.\nPASS: stable cases reach true dcgain within 1e-8.\nPASS: stable at K=47.99; unstable at K=48.01; K=48 poles -6, +/-j*sqrt(8).\nK=48 and K=50 have no convergent step final value; transient metrics are NaN.\n',residual,modalmax,maxdiff); fclose(vf);
  disp('All exp05 validation checks passed.'); diary off;
  if make_plots, plot_exp05(); end
end

function m=measure(t,y,yf)
  r=y/yf; tr=crossing(t,r,.9)-crossing(t,r,.1);
  j=find(abs(r-1)>.02,1,'last'); assert(j<numel(t));
  bnd=1+.02*sign(r(j)-1); ts=t(j)+(bnd-r(j))*(t(j+1)-t(j))/(r(j+1)-r(j));
  [peak,j]=max(r); assert(j>1 && j<numel(t));
  a=(r(j-1)-2*r(j)+r(j+1))/2; b=(r(j+1)-r(j-1))/2; d=-b/(2*a);
  peak=peak+b*d+a*d*d; tp=t(j)+d*(t(j+1)-t(j)); m=[100*(peak-1),tp,tr,ts];
end
function v=crossing(t,y,level)
  j=find(y>=level,1); assert(j>1); v=t(j-1)+(level-y(j-1))*(t(j)-t(j-1))/(y(j)-y(j-1));
end
function writecsv(path,header,A)
  f=fopen(path,'w'); fprintf(f,'%s\n',header); fclose(f); dlmwrite(path,A,'-append','precision','%.15g');
end
