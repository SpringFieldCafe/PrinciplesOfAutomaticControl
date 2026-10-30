function plot_exp03(which)
% Rebuild figures from the actual Octave calculation. Optional which=1..5.
  if nargin<1, which=1:5; end
  root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out=fullfile(root,'.local','exp03','results'); load(fullfile(out,'exp03_data.mat'));
  groups={1:4,5:7,[2 8 9],[2 10 11]};
  labels={{'zeta=0.2','zeta=0.5','zeta=0.7','zeta=1.0'}, {'wn=2 rad/s','wn=4 rad/s','wn=8 rad/s'}, {'Second-order baseline','Zero at -2 /s','Zero at -10 /s'}, {'Second-order baseline','Pole at -5 /s','Pole at -20 /s'}};
  titles={'A: Damping ratio (wn = 4 rad/s)','B: Natural frequency (zeta = 0.5)','C: Added zero, unity DC gain','C: Added pole, unity DC gain'};
  stems={'A_damping','B_frequency','C_zeros','C_poles','B_poles'};
  for n=which
    f=figure('visible','off','position',[100 100 1100 680]);
    set(0,'defaultaxesfontsize',12); set(f,'color','w');
    if n<=4
      ids=groups{n}; h=plot(t,Y(:,ids),'linewidth',1.7); hold on;
      xmax=[8 6 4 4](n);
      plot([0 xmax],[.98 .98],'k:'); plot([0 xmax],[1.02 1.02],'k:');
      plot([0 xmax],[1 1],'k--');
      for j=1:numel(ids)
        k=ids(j); if isfinite(M(k,2)), plot(M(k,2),1+M(k,1)/100,'o','color',get(h(j),'color'),'handlevisibility','off'); end
      end
      xlim([0 xmax]); ylim([0 max(max(Y(:,ids)))*1.12]);
      xlabel('Time (s)'); ylabel('Unit-step response (dimensionless)'); title(titles{n});
      legend(h,labels{n},'location','northeast');
      text(xmax*.98,.80,'Dashed: final value; dotted: 2% band','horizontalalignment','right','fontsize',10);
      text(xmax*.98,.74,'Circles: interpolated response peaks','horizontalalignment','right','fontsize',10);
    else
      hold on; c=lines(3);
      for j=1:3
        wn=[2 4 8](j); p=roots([1 wn wn^2]); plot(real(p),imag(p),'x','markersize',12,'linewidth',2,'color',c(j,:));
        text(real(p(1))+.15,imag(p(1)),sprintf('wn=%g rad/s',wn),'fontsize',11);
      end
      xlabel('Real part (1/s)'); ylabel('Imaginary part (1/s)'); title('B: Closed-loop poles, zeta = 0.5');
      xlim([-5 0]); ylim([-8 8]);
    end
    grid on; box on; set(gca,'fontsize',12);
    print(f,fullfile(out,[stems{n} '.png']),'-dpng','-r150');
    print(f,fullfile(out,[stems{n} '.svg']),'-dsvg'); close(f);
  end
end
