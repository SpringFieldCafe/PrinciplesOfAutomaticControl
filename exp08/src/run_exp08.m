function run_exp08(mode)
% GNU Octave/control experiment 08. Modes: nominal (default), robust.
% Paths anchored to source location. Unit step/ramp, zero initial states.
  if nargin<1, mode='nominal'; end
  pkg load control; root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out=fullfile(root,'.local','exp08','results'); if ~exist(out,'dir'),mkdir(out);end
  s=tf('s'); G=1/(s*(s+2)*(s+5));
  wc=2; phi=35; ratio=(1+sind(phi))/(1-sind(phi));
  z=wc/sqrt(ratio);p=wc*sqrt(ratio);k=wc*sqrt(wc^2+4)*sqrt(wc^2+25)*sqrt(ratio);
  Cs={tf(10),k*(s+z)/(s+p),48*(s+2)/(s+12),48*(s+2)/(s+12)*(s+.05)/(s+1/30)};
  pars=[0 10 0 0 0 0;1 k z p 0 0;2 48 2 12 0 0;3 48 2 12 .05 1/30];
  t=(0:.01:400)';w=logspace(-4,4,5001)';
  if strcmp(mode,'robust')
    rows=[];
    for id=[1 3]
      for gain=[.8 1 1.2]
        for a=[1.6 2 2.4]
          for b=[4 5 6]
            plant=gain/(s*(s+a)*(s+b)); L=Cs{id+1}*plant; T=feedback(L,1);
            poles=pole(T); stable=all(real(poles)<0);
            [gm,pm]=margin(L);kv=dcgain(Cs{id+1})*gain/(a*b);
            if stable, y=step(T,t);m=metrics(t,y,dcgain(T));else,m=nan(1,5);end
            passed=stable && pm>=50 && m(1)<=20 && m(3)<=4 && kv>=1;
            rows=[rows;id,gain,a,b,stable,pm,kv,m(1),m(3),max(real(poles)),passed];
          end
        end
      end
    end
    csvout(out,'robustness.csv','id,plant_gain,pole_a,pole_b,stable,PM_deg,Kv,overshoot_percent,settling_2pct_s,max_pole_real,all_targets',rows);
    for id=[1 3],r=rows(rows(:,1)==id,:);fprintf('id %d: %d/%d stable; %d/%d meet targets\n',id,sum(r(:,5)),size(r,1),sum(r(:,11)),size(r,1));end
    return;
  end
  assert(strcmp(mode,'nominal'));
  csvout(out,'parameters.csv','id,K,lead_zero,lead_pole,lag_zero,lag_pole',pars);
  tab=[];ps=[];bodes=w;responses=t;st=w;ramps=t;validation=[];
  for j=1:4
    C=Cs{j};L=C*G;T=feedback(L,1);S=feedback(tf(1),L);CS=C*S;
    pp=pole(T);assert(all(real(pp)<0));
    [gm,pm,wphase,wgain]=margin(L);kv=dcgain(C)/10;
    y=step(T,t);yf=dcgain(T);m=metrics(t,y,yf);
    f=abs(horner_response(T,w)); ix=find(f<=abs(yf)/sqrt(2),1);
    bw=exp(interp1(log(f(ix-1:ix)),log(w(ix-1:ix)),log(abs(yf)/sqrt(2))));
    lj=horner_response(L,w);sj=horner_response(S,w);tj=horner_response(T,w);uj=horner_response(CS,w);
    ms=max(abs(sj));mt=max(abs(tj));noise100=abs(squeeze(freqresp(CS,100)));
    tab=[tab;j-1,kv,gm,20*log10(gm),pm,wphase,wgain,bw,m,abs(1-yf),ms,mt,noise100,dcgain(C),pars(j,2),pm>=50 && m(1)<=20 && m(3)<=4 && kv>=1];
    ps=[ps;repmat(j-1,numel(pp),1),real(pp),imag(pp)];
    responses=[responses y];bodes=[bodes,20*log10(abs(lj)),unwrap(angle(lj))*180/pi];
    st=[st,abs(sj),abs(tj),abs(uj)];
    yr=lsim(T,t,t);er=t-yr;ramps=[ramps er];
    % Smaller output interval and exact DC/ramp predictions cross-check metrics.
    tfine=(0:.005:20)';fine=step(T,tfine);common=max(abs(fine(1:2:end)-y(1:2001)));
    fine_m=metrics(tfine,fine,yf);
    identity=max(abs(sj+tj-1));ramp_err=abs(er(end)-1/kv);
    assert(common<1e-8);assert(identity<1e-10);assert(ramp_err<1e-5);
    validation=[validation;j-1,common,abs(m(3)-fine_m(3)),ramp_err,identity];
  end
  csvout(out,'metrics.csv','id,Kv,GM_factor,GM_dB,PM_deg,phase_cross_rad_s,gain_cross_rad_s,BW_rad_s,overshoot_percent,peak_time_s,settling_2pct_s,rise_10_90_s,peak_value,step_error,Ms,Mt,CS_at_100,C_DC,C_highfreq,all_targets',tab);
  csvout(out,'poles.csv','id,real_per_s,imag_rad_s',ps);
  csvout(out,'step.csv','time_s,baseline,lead_A,lead_B,lead_lag',responses);
  csvout(out,'bode.csv','omega_rad_s,base_dB,base_deg,leadA_dB,leadA_deg,leadB_dB,leadB_deg,lag_dB,lag_deg',bodes);
  csvout(out,'sensitivity.csv','omega_rad_s,base_S,base_T,base_CS,leadA_S,leadA_T,leadA_CS,leadB_S,leadB_T,leadB_CS,lag_S,lag_T,lag_CS',st);
  csvout(out,'ramp_error.csv','time_s,baseline,lead_A,lead_B,lead_lag',ramps);
  csvout(out,'validation.csv','id,dt_common_error,settling_dt_difference,ramp_final_error,SplusT_error',validation);
  save(fullfile(out,'models.mat'),'Cs','G','pars','tab','validation');
  fid=fopen(fullfile(out,'environment.txt'),'w');fprintf(fid,'GNU Octave %s\n',version);v=pkg('list');for j=1:numel(v),if strcmp(v{j}.name,'control'),fprintf(fid,'control %s\n',v{j}.version);end;end;fprintf(fid,'All assertions PASS. dt=.01 s, horizon 400 s, fine dt=.005 s through 20 s.\n');fclose(fid);
  disp(tab);disp(validation);
end

function m=metrics(t,y,yf)
% Stable, positive-final-value systems: interpolated 10-90 rise and final
% 2 percent band entry. Horizon must include the complete settling tail.
  [peak,idx]=max(y);os=max(0,100*(peak-yf)/abs(yf));tp=t(idx);
  if os<1e-7, os=0;tp=NaN;peak=yf;end
  lo=crossing(t,y,.1*yf);hi=crossing(t,y,.9*yf);
  err=abs(y-yf);last=find(err>.02*abs(yf),1,'last');
  if isempty(last),ts=0;elseif last==numel(t),ts=NaN;else
    ts=t(last)+(t(last+1)-t(last))*(err(last)-.02*abs(yf))/(err(last)-err(last+1));end
  m=[os,tp,ts,hi-lo,peak];
end
function v=crossing(t,y,level)
  i=find(y>=level,1);if isempty(i),v=NaN;elseif i==1,v=t(1);else,v=t(i-1)+(t(i)-t(i-1))*(level-y(i-1))/(y(i)-y(i-1));end
end
function csvout(out,name,header,data)
  path=fullfile(out,name);fid=fopen(path,'w');assert(fid>=0);fprintf(fid,'%s\n',header);fclose(fid);dlmwrite(path,data,'-append','delimiter',',','precision','%.17g');
end

function h=horner_response(sys,w)
  [num,den]=tfdata(sys,'v');h=polyval(num,1i*w)./polyval(den,1i*w);
end
