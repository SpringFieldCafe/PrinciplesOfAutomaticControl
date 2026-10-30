function run_exp07(mode)
% Octave 11 + control 4.2.1; add source folder to path, then run_exp07.
% Optional 'compute' / 'plot' separates numerical work and rendering.
% Unity negative feedback, zero initial state, unit reference step.
% Signals normalized; time s, frequency rad/s. Output is project-private.
  if nargin<1,mode='all';end
  root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out=fullfile(root,'.local','exp07','results'); if ~exist(out,'dir'),mkdir(out);end
  pkg load control;
  if any(strcmp(mode,{'all','compute'}))
    Ks=[20 80 90 100]; den=[1 7 14 8]; w=unique([logspace(-4,4,12001) sqrt(14)])';
    margins=[]; poles=[]; winding=[]; steps=[]; convergence=[]; freq=[]; nyq=[]; verification=[];
    for j=1:4
      K=Ks(j); L=tf(K,den); T=feedback(L,1); p=pole(T);
      poles=[poles; K*ones(3,1) real(p) imag(p)];
      stable=all(real(p)<-1e-9); Z=sum(real(p)>1e-9); boundary=sum(abs(real(p))<=1e-9);
      [gm,pm,wpc,wgc]=margin(L);
      % All positive gain crossings: polynomial in x=w^2; no scan omissions.
      xr=roots([1 21 84 64-K*K]); x=xr(abs(imag(xr))<1e-9 & real(xr)>0);
      assert(numel(x)==1); gc=sqrt(real(x)); pc=sqrt(14);
      pma=180-sum(atan(gc./[1 2 4]))*180/pi; gma=90/K;
      if abs(pma)<1e-10,pma=0;end
      bw=NaN;
      if stable
        xb=roots([1 21 84-14*K -(8+K)^2]); xb=real(xb(abs(imag(xb))<1e-9&real(xb)>0)); assert(numel(xb)==1);bw=sqrt(xb);
      end
      margins=[margins; K gm 20*log10(gm) pm wpc wgc gma pma pc gc bw stable];
      lv=K./polyval(den,1i*w);tv=K./polyval([1 7 14 8+K],1i*w);
      for wt=[.1 1 gc pc 10 100]
        assert(abs(squeeze(freqresp(L,wt))-K/polyval(den,1i*wt))<1e-10);
      end
      freq=[freq; K*ones(size(w)) w 20*log10(abs(lv)) -sum(atan(w./[1 2 4]),2)*180/pi abs(tv)];
      % Standard RHP contour: imaginary axis -j infinity -> +j infinity.
      % N is CLOCKWISE winding count around -1, so Z=P+N.
      wf=[-flipud(w);0;w]; lf=K./polyval(den,1i*wf); z=1+lf;
      if boundary>0,N=NaN;raw=NaN; else raw=-sum(angle(z([2:end 1])./z))/(2*pi);N=round(raw);assert(N==Z);end
      [~,mi]=min(abs(1+lv)); lo=w(max(1,mi-1));hi=w(min(numel(w),mi+1));
      [wmind,dmin]=fminbnd(@(v) abs(1+K/polyval(den,1i*v)),lo,hi,optimset('TolX',1e-12));
      if K==90,wmind=pc;dmin=0;end
      winding=[winding;K 0 N Z boundary raw dmin wmind -K/90 pc];
      nyq=[nyq; K*ones(size(wf)) wf real(lf) imag(lf)];
      if stable
        tend=max(20,14/min(-real(p))); tt=(0:.002:ceil(tend))'; yy=step(T,tt); final=dcgain(T);
        m=metrics(tt,yy,final); th=(0:.001:ceil(tend))'; yh=step(T,th); mh=metrics(th,yh,final);
        steps=[steps;K final m];convergence=[convergence;K m(3:5) mh(3:5) abs(m(3:5)-mh(3:5))];
        % Independent rational residues reproduce the unit-step response.
        [rr,pp]=residue(K,conv([1 7 14 8+K],[1 0])); ya=real(exp(tt*pp.')*rr);
        analytic_err=max(abs(yy-ya)); assert(analytic_err<1e-8);
        assert(abs(yy(end)-final)/final<2e-6);
      else
        tt=(0:.002:30)';yy=step(T,tt);analytic_err=NaN;
      end
      writecsv(out,sprintf('step_K%d',K),'t_s,y',[tt yy]);
      verification=[verification;K abs(gm-gma) abs(pm-pma) abs(wpc-pc) abs(wgc-gc) analytic_err];
      assert(max(verification(end,2:5))<1e-8);
    end
    writecsv(out,'margins_all_crossings','K,GM_ratio,GM_dB,PM_deg,w_phase_cross_rad_s,w_gain_cross_rad_s,all_GM_ratio,all_PM_deg,all_phase_cross_rad_s,all_gain_cross_rad_s,stable_bandwidth_rad_s,stable',margins);
    writecsv(out,'closed_loop_poles','K,pole_real_per_s,pole_imag_rad_s',poles);
    writecsv(out,'nyquist_counts','K,P,N_clockwise,Z_RHP,imag_axis_poles,unrounded_N,minimum_distance_to_minus1,minimum_distance_frequency_rad_s,real_axis_crossing,phase_cross_frequency_rad_s',winding);
    writecsv(out,'frequency_responses','K,w_rad_s,open_magnitude_dB,open_phase_deg,closed_magnitude',freq);
    writecsv(out,'nyquist_coordinates','K,w_rad_s,real_L,imag_L',nyq);
    writecsv(out,'stable_step_metrics','K,final_value,peak,peak_time_s,rise10_90_s,settle2pct_s,overshoot_pct',steps);
    writecsv(out,'step_grid_convergence','K,rise_dt002,settle_dt002,overshoot_dt002,rise_dt001,settle_dt001,overshoot_dt001,rise_difference_s,settle_difference_s,overshoot_difference_pct',convergence);
    writecsv(out,'verification','K,GM_difference,PM_difference_deg,phase_cross_difference,gain_cross_difference,step_residue_max_error',verification);
    % Bode readout: coarse log-frequency points with interpolation, independent
    % of the polynomial solution. Save bracketing coordinates for a hand check.
    wr=logspace(-1,1,81)';mag=20*log10(abs(80./polyval(den,1i*wr)));ph=-sum(atan(wr./[1 2 4]),2)*180/pi;
    ig=find(mag<=0,1);ip=find(ph<=-180,1);
    rg=10^interp1(mag(ig-1:ig),log10(wr(ig-1:ig)),0);
    rp=10^interp1(ph(ip-1:ip),log10(wr(ip-1:ip)),-180);
    readpm=180+interp1(log10(wr),ph,log10(rg));readgm=-interp1(log10(wr),mag,log10(rp));
    writecsv(out,'bode_readout_K80','w_rad_s,magnitude_dB,phase_deg',[wr mag ph]);
    writecsv(out,'bode_readout_comparison','K,read_GM_dB,exact_GM_dB,read_PM_deg,exact_PM_deg,read_w_phase,exact_w_phase,read_w_gain,exact_w_gain',[80 readgm margins(2,3) readpm margins(2,4) rp margins(2,5) rg margins(2,6)]);
    save(fullfile(out,'workspace.mat'),'Ks','w','margins','poles','winding','steps','freq','nyq');
    fid=fopen(fullfile(out,'environment.txt'),'w');fprintf(fid,'GNU Octave %s\n',version);ps=pkg('list');for k=1:numel(ps),if strcmp(ps{k}.name,'control'),fprintf(fid,'control %s\n',ps{k}.version);end;end;fclose(fid);
    disp(margins);disp(winding);disp(steps);disp(verification);
  end
  if any(strcmp(mode,{'all','plot'}))
    load(fullfile(out,'workspace.mat'));
    set(0,'defaultfigurevisible','off','defaultaxesfontsize',11,'defaultlinelinewidth',1.6,'defaultaxescolororder',[0 .32 .65;.85 .23 .10;.20 .55 .25;.55 .25 .65]);
    labels={'K=20','K=80','K=90','K=100'};colors=get(0,'defaultaxescolororder');
    f=figure('position',[50 50 1200 650]);
    for panel=1:2
      subplot(1,2,panel);hold on;
      for j=1:4
        ix=nyq(:,1)==Ks(j);v=nyq(ix,3)+1i*nyq(ix,4);plot(real(v),imag(v),'color',colors(j,:));
      end
      plot(-1,0,'kp','markersize',12,'markerfacecolor','k');grid on;axis equal;xlabel('Real L(jw)');ylabel('Imaginary L(jw)');
      if panel==1,title('Full Nyquist locus (-infinity < w < infinity)');legend([labels {'-1'}],'location','northeast');
      else,xlim([-1.35 -.65]);ylim([-.45 .45]);title('Critical point and clockwise crossings');
        for j=1:4
          K=Ks(j);wa=sqrt(14)-.05;wb=sqrt(14)+.05;va=K/polyval([1 7 14 8],1i*wa);vb=K/polyval([1 7 14 8],1i*wb);
          quiver(real(va),imag(va),real(vb-va),imag(vb-va),0,'color',colors(j,:),'linewidth',2,'maxheadsize',1);
        end
        text(-1.33,.38,'Arrows: increasing positive w','fontsize',10);
      end
    end
    exportfig(f,out,'nyquist_full_and_local');
    f=figure('position',[50 50 1100 750]);
    hb=[];
    for j=1:4
      subplot(2,2,j);ix=nyq(:,1)==Ks(j);v=nyq(ix,3)+1i*nyq(ix,4);plot(real(v),imag(v));hold on;plot(-1,0,'rp','markersize',12,'markerfacecolor','r');
      plot(winding(j,9),0,'ko','markersize',6);grid on;xlim([-1.6 .5]);ylim([-1.1 1.1]);xlabel('Real L');ylabel('Imaginary L');
      if Ks(j)==90,title('K=90: through -1, N undefined; Z=0, 2 axis poles');else,title(sprintf('K=%d: P=0, clockwise N=%d, Z=%d',Ks(j),winding(j,3),winding(j,4)));end
    end
    exportfig(f,out,'nyquist_counts');
    f=figure('position',[50 50 1100 780]);
    for j=1:4
      ix=freq(:,1)==Ks(j);a=freq(ix,:);
      subplot(2,1,1);hb(j)=semilogx(a(:,2),a(:,3),'color',colors(j,:));hold on;
      plot(margins(j,6),0,'o','color',colors(j,:));
      subplot(2,1,2);semilogx(a(:,2),a(:,4),'color',colors(j,:));hold on;
      plot(margins(j,6),-180+margins(j,4),'o','color',colors(j,:));
    end
    subplot(2,1,1);plot([.1 100],[0 0],'k--');xlim([.1 100]);ylim([-90 25]);grid on;ylabel('Magnitude (dB)');title('Bode: circles indicate gain crossings');legend(hb,labels,'location','southwest');
    subplot(2,1,2);plot([.1 100],[-180 -180],'k--');xlim([.1 100]);grid on;ylabel('Phase (deg)');xlabel('Frequency (rad/s)');title('All four phase curves coincide');
    exportfig(f,out,'bode_comparison');
    f=figure('position',[50 50 1000 700]);a=dlmread(fullfile(out,'bode_readout_K80.csv'),',',1,0);
    subplot(2,1,1);semilogx(a(:,1),a(:,2),'o-');hold on;plot([.1 10],[0 0],'k--');plot(margins(2,5),-margins(2,3),'rp','markersize',12);xlim([2.8 4.3]);ylim([-3 5]);grid on;ylabel('Magnitude (dB)');title(sprintf('K=80: GM = %.4f dB at %.4f rad/s',margins(2,3),margins(2,5)));
    subplot(2,1,2);semilogx(a(:,1),a(:,3),'o-');hold on;plot([.1 10],[-180 -180],'k--');plot(margins(2,6),margins(2,4)-180,'rp','markersize',12);xlim([2.8 4.3]);ylim([-185 -162]);grid on;ylabel('Phase (deg)');xlabel('Frequency (rad/s)');title(sprintf('PM = %.4f deg at %.4f rad/s; circles = coarse readout grid',margins(2,4),margins(2,6)));
    for pp=1:2,subplot(2,1,pp);set(gca,'xtick',[2.8 3 3.2 3.4 3.6 3.8 4 4.2],'xticklabel',{'2.8','3.0','3.2','3.4','3.6','3.8','4.0','4.2'});end
    exportfig(f,out,'bode_K80_readout');
    f=figure('position',[50 50 1150 800]);
    for j=1:4
      a=dlmread(fullfile(out,sprintf('step_K%d.csv',Ks(j))),',',1,0);subplot(2,2,j);plot(a(:,1),a(:,2),'color',colors(j,:));hold on;grid on;xlabel('Time (s)');ylabel('Output y');
      if j<=2,fin=steps(j,2);plot([0 a(end,1)],[fin fin],'k--');plot([0 a(end,1)],[1.02*fin 1.02*fin],'k:');plot([0 a(end,1)],[.98*fin .98*fin],'k:');title(sprintf('K=%d: OS %.2f%%; settling %.3f s',Ks(j),steps(j,7),steps(j,6)));xlim([0 min(a(end,1),1.25*steps(j,6))]);
      elseif j==3,title('K=90: sustained oscillation (imaginary-axis poles)');else,title('K=100: growing oscillation (2 RHP poles)');end
    end
    exportfig(f,out,'step_all_gains');
    f=figure('position',[50 50 1000 700]);
    subplot(2,1,1);hold on;
    for j=1:2,a=dlmread(fullfile(out,sprintf('step_K%d.csv',Ks(j))),',',1,0);plot(a(:,1),a(:,2)/steps(j,2));end
    xlim([0 6]);grid on;ylabel('Output / final value');xlabel('Time (s)');title('Stable responses normalized by true DC gain');legend(labels(1:2));
    subplot(2,1,2);hold on;
    for j=1:2,ix=freq(:,1)==Ks(j);semilogx(freq(ix,2),20*log10(freq(ix,5)/steps(j,2)),'color',colors(j,:));plot(margins(j,11),-10*log10(2),'o','color',colors(j,:));end
    set(gca,'xscale','log');plot([.1 100],[-10*log10(2) -10*log10(2)],'k--');xlim([.1 100]);ylim([-80 35]);grid on;xlabel('Frequency (rad/s)');ylabel('Normalized |T| (dB)');title('Bandwidth: first downward crossing of DC / sqrt(2)');
    exportfig(f,out,'bandwidth_and_speed');
  end
end
function m=metrics(t,y,fin)
  [pk,k]=max(y);i=find(abs(y-fin)>.02*abs(fin),1,'last');if i==numel(t),ts=NaN;else,ts=t(i+1);end
  m=[pk t(k) cross(t,y,.9*fin)-cross(t,y,.1*fin) ts 100*(pk/fin-1)];
end
function tc=cross(t,y,v)
  k=find(y>=v,1);tc=t(k-1)+(v-y(k-1))*(t(k)-t(k-1))/(y(k)-y(k-1));
end
function writecsv(out,name,header,A)
  p=fullfile(out,[name '.csv']);fid=fopen(p,'w');assert(fid>=0);fprintf(fid,'%s\n',header);fclose(fid);dlmwrite(p,A,'-append','precision','%.15g');
end
function exportfig(f,out,name)
  print(f,fullfile(out,[name '.png']),'-dpng','-r150');print(f,fullfile(out,[name '.svg']),'-dsvg');close(f);
end
