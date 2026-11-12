function run_exp04()
% EXP04 stability of q(s)=s^3+3s^2+2s+K. Octave + control.
% Run after adding this source directory to the Octave path.
% Time is in seconds; unit step, zero initial state unless stated otherwise.
% All generated files remain in the project-private output directory.
  pkg load control;
  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out = fullfile(root,'.local','exp04','results');
  if ~exist(out,'dir'), mkdir(out); end
  K = [1 5.9 6 8]; poles = []; checks = [];
  t = (0:0.02:60)'; ys = zeros(numel(t),4);
  tl = (0:0.05:1200)'; yl = zeros(numel(tl),3);
  for j=1:4
    k=K(j); q=[1 3 2 k]; p=roots(q); T=tf(k,q); pp=pole(T);
    match=max(arrayfun(@(z) min(abs(z-pp)),p));
    residual=max(abs(polyval(q,p))./(abs(p).^3+3*abs(p).^2+2*abs(p)+abs(k)));
    tol=1e-10*max(1,max(abs(p))); nr=sum(real(p)>tol);
    first=[1 3 (6-k)/3 k];
    if k==6, changes=0; else, changes=sum(first(1:end-1).*first(2:end)<0); end
    poles=[poles; repmat(k,3,1),real(p),imag(p)];
    checks=[checks;k,max(real(p)),nr,changes,match,residual,tol];
    assert(nr==changes); assert(match<1e-10); assert(residual<1e-12);
    ys(:,j)=step(T,t);
    if j<=3, yl(:,j)=step(T,tl); end
  end
  csvout(out,'poles.csv','K,real_per_s,imag_rad_per_s',poles);
  csvout(out,'checks.csv','K,max_real,RHP_roots,Routh_changes,roots_pole_error,normalized_residual,axis_tolerance',checks);
  csvout(out,'step_overlay.csv','time_s,K1,K5p9,K6,K8',[t ys]);
  csvout(out,'step_long.csv','time_s,K1,K5p9,K6',[tl yl]);
  % Repeated scans include the exact theoretical boundary, with raw and
  % tolerance-aware classifications retained separately.
  hs=[0.05 0.005 0.0005]; summary=[];
  for n=1:numel(hs)
    h=hs(n); ks=linspace(4,8,round(4/h)+1)'; a=zeros(size(ks));
    tracks=[];
    for j=1:numel(ks)
      p=roots([1 3 2 ks(j)]); a(j)=max(real(p));
      if n==1, tracks=[tracks;repmat(ks(j),3,1),real(p),imag(p)]; end
    end
    raw=find(a>0,1); neg=find(a < -1e-10,1,'last'); pos=find(a>1e-10,1);
    lo=ks(neg); hi=ks(pos); estimate=(lo+hi)/2;
    summary=[summary;h,ks(raw-1),ks(raw),lo,hi,estimate,abs(estimate-6),(hi-lo)/2];
    csvout(out,sprintf('scan_%d.csv',n),'K,max_real,raw_sign,tolerance_sign',[ks,a,sign(a),(a>1e-10)-(a < -1e-10)]);
    if n==1, csvout(out,'pole_trajectories.csv','K,real_per_s,imag_rad_per_s',tracks); end
  end
  csvout(out,'scan_summary.csv','step,raw_left,raw_first_positive,tol_left,tol_right,estimate,true_error,bracket_halfwidth',summary);
  % Independent shifted bracket avoids accidental exact-grid success.
  lo=5.973; hi=6.021; history=[];
  for j=1:35
    mid=(lo+hi)/2; a=max(real(roots([1 3 2 mid])));
    history=[history;j,lo,hi,mid,a];
    if a>0, hi=mid; else, lo=mid; end
  end
  csvout(out,'bisection.csv','iteration,left,right,midpoint,max_real',history);
  assert(abs((lo+hi)/2-6)<1e-10);
  % Exact K=6 solution verifies dynamics and persistent oscillation.
  exact=1-2/11*exp(-3*t)-9/11*cos(sqrt(2)*t)-3*sqrt(2)/11*sin(sqrt(2)*t);
  analytic_error=max(abs(ys(:,3)-exact));
  th=(0:0.01:60)'; fine=step(tf(6,[1 3 2 6]),th);
  sampling_error=max(abs(fine(1:2:end)-ys(:,3)));
  assert(analytic_error<1e-9); assert(sampling_error<1e-9);
  csvout(out,'critical_exact.csv','time_s,numeric,analytic,error',[t,ys(:,3),exact,ys(:,3)-exact]);
  % K=0: x=[y ydot yddot], B=0, C=[1 0 0].
  % A has a zero eigenvalue, invisible in the zero-state transfer T=0.
  A=[0 1 0;0 0 1;0 -2 -3]; B=zeros(3,1); C=[1 0 0];
  system0=ss(A,B,C,0); ti=(0:0.05:20)';
  zi=initial(system0,[1;0;0],ti); zs=step(system0,ti);
  assert(max(abs(zi-1))<1e-10); assert(max(abs(zs))<1e-12);
  csvout(out,'K0_internal.csv','time_s,unit_initial_y,zero_state_step',[ti,zi,zs]);
  offsets=[-1e-8 -1e-10 -1e-12 0 1e-12 1e-10 1e-8]; sensitivity=[];
  for d=offsets
    p=roots([1 3 2 6+d]); a=max(real(p));
    sensitivity=[sensitivity;d,a,d/22,(a>1e-10)-(a < -1e-10)];
  end
  csvout(out,'roundoff.csv','K_offset,max_real,linear_prediction,tolerance_sign',sensitivity);
  fid=fopen(fullfile(out,'validation.txt'),'w');
  fprintf(fid,'GNU Octave %s; control package loaded\n',version);
  v=pkg('list'); for n=1:numel(v), if strcmp(v{n}.name,'control'), fprintf(fid,'control %s\n',v{n}.version); end; end
  fprintf(fid,'Unit step, zero initial conditions; time s; K dimensionless.\n');
  fprintf(fid,'analytic K6 max error %.17g\nsampling dt 0.02 vs 0.01 max error %.17g\n',analytic_error,sampling_error);
  fprintf(fid,'bisection final [%.17g, %.17g]; midpoint error %.17g\n',lo,hi,abs((lo+hi)/2-6));
  fprintf(fid,'All assertions PASS.\n'); fclose(fid);
  save(fullfile(out,'exp04_data.mat'),'K','poles','checks','t','ys','tl','yl','summary','history','sensitivity');
  disp(checks); fprintf('PASS: exp04 numerical evidence saved.\n');
end

function csvout(out,name,header,data)
  path=fullfile(out,name); fid=fopen(path,'w'); assert(fid>=0);
  fprintf(fid,'%s\n',header); fclose(fid);
  dlmwrite(path,data,'-append','delimiter',',','precision','%.17g');
end
