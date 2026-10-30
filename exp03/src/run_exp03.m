function run_exp03(make_plots)
% GNU Octave 11 + control. Zero initial state; unit step; seconds, rad/s.
% Entry: addpath(fullfile(project_root,'exp03','src')); run_exp03
% Outputs are private and resolved from this source, independently of pwd.
  if nargin < 1, make_plots = true; end
  pkg load control;
  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out = fullfile(root,'.local','exp03','results');
  if ~exist(out,'dir'), mkdir(out); end
  diary(fullfile(out,'run.log')); diary on;
  printf('GNU Octave %s; zero initial state; unit step; dt=0.0005 s; horizon=30 s\n',version);
  pkg list control;
  z=[.2 .5 .7 1 .5 .5 .5 .5 .5 .5 .5];
  w=[4 4 4 4 2 4 8 4 4 4 4];
  names={'A_z02','A_z05','A_z07','A_z10','B_w2','B_w4','B_w8','C_zero2','C_zero10','C_pole5','C_pole20'};
  t=(0:.0005:30)'; tc=t(1:2:end); Y=zeros(numel(t),11);
  fid=fopen(fullfile(out,'metrics.csv'),'w');
  fprintf(fid,'case,zeta,wn_rad_s,added_zero,added_pole,dcgain,Mp_percent,tp_s,tr10_90_s,ts2pct_s,Mp_theory_percent,tp_theory_s,tr_theory_s,ts_exact_theory_s,ts_approx_s,Mp_error_pp,tp_error_s,tr_error_s,ts_error_s,analytic_waveform_max_error,final_error,dt_Mp_delta_pp,dt_tp_delta_s,dt_tr_delta_s,dt_ts_delta_s\n');
  pf=fopen(fullfile(out,'poles.csv'),'w'); fprintf(pf,'case,pole_real_per_s,pole_imag_per_s\n');
  M=zeros(11,4); checks=zeros(11,4);
  for k=1:11
    sys=tf(w(k)^2,[1 2*z(k)*w(k) w(k)^2]); az=NaN; ap=NaN;
    if k==8, az=-2; elseif k==9, az=-10; elseif k==10, ap=-5; elseif k==11, ap=-20; end
    if isfinite(az), sys=sys*tf([1/(-az) 1],1); end
    if isfinite(ap), sys=sys*tf(-ap,[1 -ap]); end
    p=pole(sys); assert(all(real(p)<0));
    for j=1:numel(p), fprintf(pf,'%s,%.15g,%.15g\n',names{k},real(p(j)),imag(p(j))); end
    yf=dcgain(sys); assert(abs(yf-1)<1e-12);
    y=step(sys,t); yc=step(sys,tc); Y(:,k)=y;
    mono=(k==4); m=measure(t,y,yf,mono); mc=measure(tc,yc,yf,mono); M(k,:)=m;
    delta=abs(m-mc); if mono, delta(2)=0; end; checks(k,:)=delta;
    assert(delta(1)<.001 && delta(2)<1e-5 && delta(3)<1e-5 && delta(4)<1e-5);
    assert(abs(y(end)-yf)<1e-8);
    th=NaN(1,4); approx=NaN; wave=NaN;
    if k<=7
      f=@(x) analytic(x,z(k),w(k)); wave=max(abs(y-f(t))); assert(wave<1e-9);
      th(1)=0; th(2)=NaN;
      if z(k)<1, th(1)=100*exp(-pi*z(k)/sqrt(1-z(k)^2)); th(2)=pi/(w(k)*sqrt(1-z(k)^2)); end
      th(3)=exact_cross(f,t,.9)-exact_cross(f,t,.1);
      q=find(abs(f(t)-1)>.02,1,'last');
      th(4)=fzero(@(x) abs(f(x)-1)-.02,[t(q),t(q+1)]);
      approx=4/(z(k)*w(k));
      err=abs(m-th); if mono, err(2)=0; end
      assert(err(1)<.001 && max(err(2:4))<1e-5);
    end
    row=[z(k),w(k),az,ap,yf,m,th,approx,m-th,wave,abs(y(end)-yf),delta];
    fprintf(fid,'%s',names{k}); fprintf(fid,',%.15g',row); fprintf(fid,'\n');
    printf('%s: Mp=%.6f%% tp=%.6f tr=%.6f ts=%.6f s\n',names{k},m);
  end
  fclose(fid); fclose(pf);
  write_csv(fullfile(out,'responses.csv'),['time_s,' strjoin(names,',')],[t,Y]);
  ef=fopen(fullfile(out,'second_order_approximation.csv'),'w');
  fprintf(ef,'case,max_abs_response_error,ISE_s,Mp_delta_pp,tp_delta_s,tr_delta_s,ts_delta_s\n');
  for k=8:11
    e=Y(:,k)-Y(:,2); fprintf(ef,'%s',names{k}); fprintf(ef,',%.15g',[max(abs(e)),trapz(t,e.^2),M(k,:)-M(2,:)]); fprintf(ef,'\n');
  end
  fclose(ef);
  % Same damping gives exactly the same response in normalized time wn*t.
  tau=(0:.002:40)'; normY=zeros(numel(tau),3);
  for j=1:3, k=4+j; normY(:,j)=step(tf(w(k)^2,[1 w(k) w(k)^2]),tau/w(k)); end
  scale_error=max(max(abs(normY-normY(:,1)))); assert(scale_error<1e-10);
  write_csv(fullfile(out,'normalized_time.csv'),'wn_times_t,response_w2,response_w4,response_w8',[tau,normY]);
  save(fullfile(out,'exp03_data.mat'),'t','Y','M','names','z','w');
  vf=fopen(fullfile(out,'validation.txt'),'w');
  fprintf(vf,'PASS: stable poles, dcgain=1, final error <1e-8 for all 11 cases.\nPASS: all 7 standard second-order curves vs analytic solution <1e-9.\nPASS: standard metrics vs analytic formulas/root solving: Mp <0.001 pp; times <1e-5 s.\nPASS: dt 0.001 -> 0.0005 s: maximum metric changes %.9g pp, %.9g s, %.9g s, %.9g s.\nPASS: normalized time response error %.9g.\nCritical damping: no finite peak; tp=NaN, Mp=0.\n',max(checks,[],1),scale_error); fclose(vf);
  disp('All numerical validation checks passed.'); diary off;
  if make_plots, plot_exp03(); end
end

function y=analytic(t,z,w)
  if z==1, y=1-(1+w*t).*exp(-w*t);
  else, wd=w*sqrt(1-z*z); y=1-exp(-z*w*t).*(cos(wd*t)+z/sqrt(1-z*z)*sin(wd*t)); end
end

function v=exact_cross(f,t,level)
  j=find(f(t)>=level,1); v=fzero(@(x) f(x)-level,[t(j-1),t(j)]);
end

function m=measure(t,y,yf,monotone)
% Linear crossing interpolation; quadratic peak interpolation; true dcgain.
  r=y/yf; tr=crossing(t,r,.9)-crossing(t,r,.1);
  q=find(abs(r-1)>.02,1,'last'); assert(q<numel(t));
  if isempty(q), ts=0; else
    bound=1+.02*sign(r(q)-1); ts=t(q)+(bound-r(q))*(t(q+1)-t(q))/(r(q+1)-r(q));
  end
  if monotone, mp=0; tp=NaN;
  else
    [peak,j]=max(r); assert(j>1 && j<numel(t));
    a=(r(j-1)-2*r(j)+r(j+1))/2; b=(r(j+1)-r(j-1))/2;
    d=-b/(2*a); tp=t(j)+d*(t(j+1)-t(j)); peak=peak+b*d+a*d*d;
    mp=max(0,100*(peak-1));
  end
  m=[mp,tp,tr,ts];
end

function v=crossing(t,y,level)
  j=find(y>=level,1); assert(j>1); v=t(j-1)+(level-y(j-1))*(t(j)-t(j-1))/(y(j)-y(j-1));
end

function write_csv(path,header,A)
  f=fopen(path,'w'); fprintf(f,'%s\n',header); fclose(f); dlmwrite(path,A,'-append','precision','%.15g');
end
