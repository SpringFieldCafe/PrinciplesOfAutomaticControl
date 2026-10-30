function run_exp02(mode)
% Feedback characteristics, GNU Octave 11 + control 4.2.1.
% Entry: add this source directory to path; run_exp02('compute'), then
% run_exp02('plot'). All paths are derived from this source location.
% Time in seconds, angular frequency rad/s; signals use normalized units.
% Zero initial conditions. Negative unity feedback; d at plant input,
% n added to measured output. u = K*(r-y-n), y = G*(u+d).
  if nargin < 1, mode = 'all'; end
  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out = fullfile(root,'.local','exp02','results');
  if ~exist(out,'dir'), mkdir(out); end
  pkg load control;
  if strcmp(mode,'all') || strcmp(mode,'compute')
    s=tf('s'); G=1/(s*(s+2)); Ks=[1 4 10]; av=[1.6 1.8 2 2.2 2.4];
    w=logspace(-2,2,801)'; t=(0:0.002:20)'; tn=(0:0.001:20)';
    metrics=[]; noise=[]; freq=[]; sens=[]; check=[]; convergence=[];
    yr=[]; yd=[]; ur=[]; ud=[]; yn=[]; un=[]; sm=[]; tm=[]; dm=[];
    for j=1:3
      K=Ks(j); L=K*G; T=feedback(L,1); S=feedback(1,L); D=minreal(G*S);
      assert(all(real(pole(T))<0));
      yr(:,j)=step(T,t); yd(:,j)=step(D,t);
      ur(:,j)=step(K*S,t); ud(:,j)=step(-T,t);
      ref=measure(t,yr(:,j),dcgain(T)); dst=measure(t,yd(:,j),dcgain(D));
      th=(0:0.001:20)'; rh=measure(th,step(T,th),dcgain(T));
      convergence=[convergence; K ref(3:5) rh(3:5) abs(ref(3:5)-rh(3:5))];
      if K==1, analytic=1-(1+t).*exp(-t); else b=sqrt(K-1); analytic=1-exp(-t).*(cos(b*t)+sin(b*t)/b); end
      assert(max(abs(analytic-yr(:,j)))<1e-9);
      metrics=[metrics; K dcgain(T) ref dcgain(D) dst max(abs(ur(:,j)))];
      nv=0.1*sin(30*tn); yn(:,j)=lsim(-T,nv,tn); un(:,j)=lsim(-K*S,nv,tn);
      ix=tn>=12; B=[sin(30*tn(ix)) cos(30*tn(ix)) ones(sum(ix),1)];
      cy=B\yn(ix,j); cu=B\un(ix,j); ya=norm(cy(1:2)); ua=norm(cu(1:2));
      exacty=0.1*abs(evalfr_local(-T,30)); exactu=0.1*abs(evalfr_local(-K*S,30));
      fine=(0:0.0005:20)'; yf=lsim(-T,0.1*sin(30*fine),fine); uf=lsim(-K*S,0.1*sin(30*fine),fine);
      iz=fine>=12; Bf=[sin(30*fine(iz)) cos(30*fine(iz)) ones(sum(iz),1)];
      cf=Bf\yf(iz); cuf=Bf\uf(iz);
      noise=[noise; K exacty ya norm(cf(1:2)) exactu ua norm(cuf(1:2)) max(abs(un(:,j)))];
      sv=squeeze(freqresp(S,w)); tv=squeeze(freqresp(T,w)); dv=squeeze(freqresp(D,w));
      sm(:,j)=abs(sv); tm(:,j)=abs(tv); dm(:,j)=abs(dv);
      check=[check; K max(abs(sv+tv-1)) max(abs(yr(:,j)-K*yd(:,j))) abs(yr(end,j)-1) abs(yd(end,j)-1/K) abs(ya-exacty)/exacty];
      assert(check(end,2)<1e-10 && check(end,3)<1e-9 && check(end,4)<1e-6);
      assert(abs(norm(cf(1:2))-exacty)/exacty<3e-5);
      assert(abs(norm(cuf(1:2))-exactu)/exactu<3e-5);
      for a=av
        Ga=1/(s*(s+a)); La=K*Ga; Ta=feedback(La,1);
        lv=squeeze(freqresp(La,w)); ta=squeeze(freqresp(Ta,w));
        freq=[freq; [K*ones(size(w)) a*ones(size(w)) w 20*log10(abs(lv)) unwrap(angle(lv))*180/pi 20*log10(abs(ta)) unwrap(angle(ta))*180/pi]];
        for ww=[0.1 1 2 5 30]
          l0=evalfr_local(L,ww); t0=evalfr_local(T,ww); la=evalfr_local(La,ww); tt=evalfr_local(Ta,ww);
          sens=[sens; K a ww abs(la) angle(la)*180/pi abs(tt) angle(tt)*180/pi 100*(abs(la)/abs(l0)-1) angle(la/l0)*180/pi 100*(abs(tt)/abs(t0)-1) angle(tt/t0)*180/pi];
        end
      end
    end
    writecsv(out,'step_metrics','K,reference_final,reference_peak,reference_peak_time_s,reference_rise10_90_s,reference_settle2pct_s,reference_overshoot_pct,disturbance_final,disturbance_peak,disturbance_peak_time_s,disturbance_rise10_90_s,disturbance_settle2pct_s,disturbance_overshoot_pct,reference_control_peak',metrics);
    writecsv(out,'noise_metrics','K,y_amplitude_exact,y_amplitude_dt001,y_amplitude_dt0005,u_amplitude_exact,u_amplitude_dt001,u_amplitude_dt0005,u_peak_transient',noise);
    writecsv(out,'verification','K,SplusT_residual,TminusKD_residual,reference_endpoint_error,disturbance_endpoint_error,noise_coarse_relative_error',check);
    writecsv(out,'step_grid_convergence','K,rise_dt002,settle_dt002,overshoot_dt002,rise_dt001,settle_dt001,overshoot_dt001,rise_difference_s,settle_difference_s,overshoot_difference_pct',convergence);
    writecsv(out,'frequency_parameter_sweep','K,a,w_rad_s,open_magnitude_db,open_phase_deg,closed_magnitude_db,closed_phase_deg',freq);
    writecsv(out,'parameter_selected_frequencies','K,a,w_rad_s,open_magnitude,open_phase_deg,closed_magnitude,closed_phase_deg,open_magnitude_change_pct,open_phase_change_deg,closed_magnitude_change_pct,closed_phase_change_deg',sens);
    writecsv(out,'time_responses','t_s,r_K1,r_K4,r_K10,d_K1,d_K4,d_K10,ur_K1,ur_K4,ur_K10,ud_K1,ud_K4,ud_K10',[t yr yd ur ud]);
    writecsv(out,'noise_responses','t_s,n,y_K1,y_K4,y_K10,u_K1,u_K4,u_K10',[tn 0.1*sin(30*tn) yn un]);
    writecsv(out,'sensitivity_functions','w_rad_s,S_K1,S_K4,S_K10,T_K1,T_K4,T_K10,Td_K1,Td_K4,Td_K10',[w sm tm dm]);
    save(fullfile(out,'workspace.mat'),'Ks','av','w','t','tn','yr','yd','ur','ud','yn','un','sm','tm','dm','freq','metrics','noise');
    fid=fopen(fullfile(out,'environment.txt'),'w'); fprintf(fid,'GNU Octave %s\n',version); p=pkg('list'); for k=1:numel(p), if strcmp(p{k}.name,'control'), fprintf(fid,'control %s\n',p{k}.version); end; end; fclose(fid);
    disp(metrics); disp(noise); disp(check);
  end
  if strcmp(mode,'all') || strcmp(mode,'plot')
    load(fullfile(out,'workspace.mat'));
    set(0,'defaultfigurevisible','off','defaultaxesfontsize',11,'defaultlinelinewidth',1.8,'defaultaxescolororder',[0 .32 .65;.85 .23 .10;.20 .55 .25;.55 .25 .65;.6 .4 .1]);
    labels={'K = 1','K = 4','K = 10'};
    f=figure('position',[50 50 1100 760]);
    subplot(2,2,1); plot(t,yr); grid on; xlim([0 10]); xlabel('Time (s)'); ylabel('Output y'); title('Unit reference step'); legend(labels,'location','southeast');
    subplot(2,2,2); plot(t,yd); grid on; xlim([0 10]); xlabel('Time (s)'); ylabel('Output y'); title('Unit plant-input disturbance: final = 1/K'); legend(labels);
    subplot(2,2,3); plot(t,ur); grid on; xlim([0 10]); xlabel('Time (s)'); ylabel('Control u'); title('Control due to reference');
    subplot(2,2,4); plot(t,ud); grid on; xlim([0 10]); xlabel('Time (s)'); ylabel('Control u'); title('Control due to disturbance');
    exportfig(f,out,'step_and_control');
    f=figure('position',[50 50 1100 750]);
    subplot(2,2,1); semilogx(w,20*log10(sm)); grid on; xlabel('Frequency (rad/s)'); ylabel('|S| (dB)'); title('Sensitivity'); legend(labels);
    subplot(2,2,2); semilogx(w,20*log10(tm)); grid on; xlabel('Frequency (rad/s)'); ylabel('|T| (dB)'); title('Complementary sensitivity');
    subplot(2,2,3); semilogx(w,20*log10(dm)); grid on; xlabel('Frequency (rad/s)'); ylabel('|Td| (dB)'); title('Disturbance channel G S');
    subplot(2,2,4); semilogx(w,sm); hold on; plot([w(1) w(end)],[1 1],'k--'); grid on; xlabel('Frequency (rad/s)'); ylabel('|S|'); title('Parameter sensitivity reduced where |S| < 1');
    exportfig(f,out,'sensitivity_functions');
    f=figure('position',[50 50 1100 750]);
    subplot(2,2,1); plot(tn,yn); grid on; xlim([0 5]); xlabel('Time (s)'); ylabel('Output noise'); title('Noise startup (zero initial state)'); legend(labels);
    subplot(2,2,2); plot(tn,un); grid on; xlim([0 1]); xlabel('Time (s)'); ylabel('Control noise'); title('Noise startup control');
    subplot(2,2,3); plot(tn,yn); grid on; xlim([18 18.6]); ylim([-0.0013 0.0013]); xlabel('Time (s)'); ylabel('Output noise'); title('Steady output: n = 0.1 sin(30t)');
    subplot(2,2,4); plot(tn,un); grid on; xlim([18 18.6]); xlabel('Time (s)'); ylabel('Control noise'); title('Steady control: transfer = -K S');
    exportfig(f,out,'measurement_noise');
    for j=1:3
      f=figure('position',[50 50 1100 750]);
      for z=1:4
        subplot(2,2,z); hold on;
        for a=av, ix=freq(:,1)==Ks(j)&abs(freq(:,2)-a)<1e-9; semilogx(freq(ix,3),freq(ix,z+3)); end
        set(gca,'xscale','log'); grid on; xlabel('Frequency (rad/s)');
        if mod(z,2)==1, ylabel('Magnitude (dB)'); else ylabel('Phase (deg)'); end
        titles={'Open loop K G','Open loop K G','Closed loop T','Closed loop T'}; title(sprintf('%s, K = %g',titles{z},Ks(j)));
        legend({'a=1.6','a=1.8','a=2.0','a=2.2','a=2.4'},'location','southwest');
      end
      exportfig(f,out,sprintf('parameter_sweep_K%d',Ks(j)));
    end
  end
end

function v=evalfr_local(sys,w)
  v=squeeze(freqresp(sys,w));
end

function m=measure(t,y,final)
% Peak sampled on 2 ms grid. Rise uses linearly interpolated first crossings.
% Settling = first sample after last violation of +/-2% of true dcgain.
% Stable systems only; simulation endpoint must already lie inside the band.
  [peak,k]=max(y); tp=t(k);
  if peak<final && all(diff(y)>=-1e-12), peak=final; tp=Inf; end
  i=find(abs(y-final)>0.02*abs(final),1,'last');
  if isempty(i), ts=0; elseif i==numel(t), ts=NaN; else ts=t(i+1); end
  m=[peak tp crossing(t,y,0.9*final)-crossing(t,y,0.1*final) ts max(0,100*(peak/final-1))];
end

function tc=crossing(t,y,level)
  k=find(y>=level,1); if isempty(k), tc=NaN; elseif k==1, tc=t(1); else tc=t(k-1)+(level-y(k-1))*(t(k)-t(k-1))/(y(k)-y(k-1)); end
end

function writecsv(out,name,header,A)
  path=fullfile(out,[name '.csv']); fid=fopen(path,'w'); assert(fid>=0); fprintf(fid,'%s\n',header); fclose(fid); dlmwrite(path,A,'-append','precision','%.15g');
end

function exportfig(f,out,name)
  print(f,fullfile(out,[name '.png']),'-dpng','-r150');
  print(f,fullfile(out,[name '.svg']),'-dsvg'); close(f);
end
