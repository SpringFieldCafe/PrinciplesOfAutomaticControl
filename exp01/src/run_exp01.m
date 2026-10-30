function run_exp01()
  % Experiment 01: mass-spring-damper models and unity negative feedback.
  % GNU Octave 11.3 + control 4.2.1. SI units: kg, N*s/m, N/m, s, m.
  % Entry: add this source directory to path, then run_exp01.
  % All output paths are resolved from this file, independent of cwd.
  pkg load control;
  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out = fullfile(root, '.local', 'exp01', 'results');
  if ~exist(out, 'dir'), mkdir(out); end
  diary(fullfile(out, 'octave_run.txt'));
  unwind_protect
    fprintf('Octave %s\n', version()); pkg list control;
    m=1; b=2; k=5; dt=0.001; t=(0:dt:10)';
    A=[0,1;-k/m,-b/m]; B=[0;1/m]; C=[1,0]; D=0;
    G=tf(1,[m,b,k]); Gz=zpk([],roots([m,b,k]),1/m);
    Gss=ss(A,B,C,D); Gauto=ss(G);
    models={G,Gz,Gss,tf(Gss),ss(Gz),tf(Gauto),zpk(G),zpk(Gss),tf(Gz)};
    names={'tf','zpk','physical_ss','ss_to_tf','zpk_to_ss','tf_ss_tf','tf_to_zpk','ss_to_zpk','zpk_to_tf'};
    Y=zeros(length(t),length(models)); eq=zeros(length(models),4);
    reference=step(G,t); pref=sort(pole(G));
    for j=1:length(models)
      Y(:,j)=step(models{j},t);
      eq(j,:)=[j,dcgain(models{j}),max(abs(sort(pole(models{j}))-pref)),max(abs(Y(:,j)-reference))];
      fprintf('%s: dcgain %.15g, pole error %.3g, step error %.3g\n',names{j},eq(j,2:4));
    end
    assert(max(eq(:,3))<1e-10); assert(max(eq(:,4))<1e-10);
    assert(max(abs(eq(:,2)-0.2))<1e-12);
    csvwrite_header(fullfile(out,'model_equivalence.csv'),'model_id,dcgain_m_per_N,pole_max_error,step_max_error_m',eq);
    csvwrite_header(fullfile(out,'open_loop_steps.csv'),['time_s,' strjoin(names,',')],[t,Y]);
    dlmwrite(fullfile(out,'state_A.csv'),A,',','precision','%.17g');
    dlmwrite(fullfile(out,'state_B.csv'),B,',','precision','%.17g');
    dlmwrite(fullfile(out,'state_C.csv'),C,',','precision','%.17g');
    [Ac,Bc,Cc,Dc]=ssdata(Gauto);
    fid=fopen(fullfile(out,'model_parameters.txt'),'w');
    fprintf(fid,'m=1 kg; b=2 N*s/m; k=5 N/m; x(0)=xdot(0)=0.\n');
    fprintf(fid,'m*xddot+b*xdot+k*x=u. G(s)=1/(s^2+2s+5).\n');
    fprintf(fid,'zpk: no finite zeros; poles -1 +/- 2j; gain=1.\n');
    fprintf(fid,'Physical state: [x;xdot], A=[0,1;-5,-2], B=[0;1], C=[1,0], D=0.\n');
    fprintf(fid,'Model IDs: %s\n',strjoin(names,', '));
    fprintf(fid,'Octave control implements zpk construction using tf/ss storage, not a distinct zpk class.\n');
    fclose(fid);
    dlmwrite(fullfile(out,'converted_state_A.csv'),Ac,',','precision','%.17g');
    dlmwrite(fullfile(out,'converted_state_B.csv'),Bc,',','precision','%.17g');
    dlmwrite(fullfile(out,'converted_state_C.csv'),Cc,',','precision','%.17g');
    [wn,zeta,p]=damp(G); disp(G); disp(Gz); disp(Gss); disp(Gauto);
    fprintf('Open-loop damp results [wn,zeta,real(p),imag(p)]\n'); disp([wn,zeta,real(p),imag(p)]);
    gains=[1,4,10]; closed=zeros(length(t),3); metrics=zeros(3,12);
    for j=1:3
      K=gains(j); L=series(tf(K),G); T=feedback(L,1);
      % Independent interconnections: parallel branches and sensitivity identity.
      Lparallel=parallel(series(tf(K/2),G),series(tf(K/2),G));
      assert(norm(freqresp(minreal(Lparallel-L),[0.1,1,10])(:))<1e-10);
      S=feedback(tf(1),L);
      assert(norm(freqresp(minreal(parallel(T,S)-tf(1)),[0.1,1,10])(:))<1e-10);
      p=sort(pole(T)); closed(:,j)=step(T,t); yinf=K/(k+K);
      wd=sqrt(k+K-1); analytic=yinf*(1-exp(-t).*(cos(wd*t)+sin(wd*t)/wd));
      err=max(abs(closed(:,j)-analytic)); assert(err<1e-10);
      [peak,ii]=max(closed(:,j)); mp=100*(peak/yinf-1);
      t10=crossing(t,closed(:,j),0.1*yinf); t90=crossing(t,closed(:,j),0.9*yinf);
      last=find(abs(closed(:,j)-yinf)>0.02*yinf,1,'last');
      if isempty(last), ts=0; elseif last==numel(t), ts=NaN; else
        ts=t(last)+(t(last+1)-t(last))*(abs(closed(last,j)-yinf)-0.02*yinf)/(abs(closed(last,j)-yinf)-abs(closed(last+1,j)-yinf));
      end
      metrics(j,:)=[K,real(p(1)),imag(p(1)),real(p(2)),imag(p(2)),yinf,mp,t(ii),t90-t10,ts,err,1/sqrt(k+K)];
    end
    csvwrite_header(fullfile(out,'closed_loop_steps.csv'),'time_s,K1,K4,K10',[t,closed]);
    csvwrite_header(fullfile(out,'closed_loop_metrics.csv'),'K,p1_real,p1_imag,p2_real,p2_imag,steady_m,overshoot_pct,peak_time_s,rise_10_90_s,settling_2pct_s,analytic_max_error_m,zeta',metrics);
    fprintf('K / steady value / overshoot pct / peak s / rise s / settling s\n'); disp(metrics(:,[1,6:10]));
    fprintf('PASS: model conversion, series/parallel/feedback and analytic response checks.\n');
  unwind_protect_cleanup
    diary off;
  end_unwind_protect
end

function tc=crossing(t,y,level)
  i=find(y>=level,1); if isempty(i), tc=NaN; elseif i==1, tc=t(1); else
    tc=t(i-1)+(t(i)-t(i-1))*(level-y(i-1))/(y(i)-y(i-1));
  end
end

function csvwrite_header(path,header,data)
  f=fopen(path,'w'); assert(f>=0); fprintf(f,'%s\n',header); fclose(f);
  dlmwrite(path,data,'-append','delimiter',',','precision','%.17g');
end
