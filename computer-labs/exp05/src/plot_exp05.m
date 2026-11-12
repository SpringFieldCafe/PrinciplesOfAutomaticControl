function plot_exp05(which)
% Figures use computed data; theoretical sketch is independently constructed
% from real-axis and complex-branch equations, clearly identified as such.
  if nargin<1, which=1:5; end
  pkg load control;
  root=fileparts(fileparts(fileparts(mfilename('fullpath')))); out=fullfile(root,'.local','exp05','results');
  load(fullfile(out,'exp05_data.mat'));
  stems={'A_theoretical_sketch','B_rlocus','C_stable_steps','C_boundary_and_unstable','B_selected_prediction'};
  for n=which
    f=figure('visible','off','position',[100 100 1100 700]); set(f,'color','w'); hold on;
    if n<=2
      if n==1
        % Im(q)=0 gives y^2=3*x^2+12*x+8, positive-gain curved branches.
        xb=-2+2/sqrt(3); x=linspace(xb,1.5,600); yy=sqrt(max(0,3*x.^2+12*x+8));
        plot(x,yy,'b','linewidth',2); plot(x,-yy,'b','linewidth',2);
        plot([-9 -4],[0 0],'b','linewidth',2); plot([-2 0],[0 0],'b','linewidth',2);
        quiver(-5,0,-.8,0,0,'b','linewidth',1.5); quiver(-.2,0,-.35,0,0,'b','linewidth',1.5);
        quiver(-1.8,0,.35,0,0,'b','linewidth',1.5);
        title('A: Theory-derived sketch (program typeset, K >= 0)');
      else
        rlocus(tf(1,[1 6 8 0])); hold on;
        title('B: Octave rlocus with zeta = 0.5 grid');
      end
      xlim([-9 2]); ylim([-6 6]);
      sgrid(.5,[]); hold on;
      for theta=[pi/3 pi 5*pi/3]
        u=linspace(0,10,100); plot(-2+u*cos(theta),u*sin(theta),'k--','linewidth',1);
      end
      plot([-4 -2 0],[0 0 0],'kx','markersize',11,'linewidth',2);
      plot(-2,0,'ks','markersize',7,'linewidth',1.5);
      plot(-2+2/sqrt(3),0,'mo','markersize',8,'linewidth',1.5);
      plot([0 0],[-sqrt(8) sqrt(8)],'ro','markersize',9,'linewidth',2);
      plot(real(selected),imag(selected),'go','markersize',9,'linewidth',2);
      text(-8.6,5.3,'Dashed asymptotes: 60, 180, 300 deg');
      text(-8.6,4.65,'Centroid: -2; real locus: (-Inf,-4), (-2,0)');
      text(-8.6,4.0,'Breakaway: -0.84530, K=3.07920');
      text(.15,3.2,'K=48'); text(-.5,1.2,sprintf('K=%.5f',Ksel));
      if n==2, legend('off'); end
      text(-8.6,-4.7,'x: open-loop poles; red circles: imaginary-axis crossing');
      text(-8.6,-5.3,'Green circle: selected point; dotted rays: zeta=0.5');
      xlabel('Real part (1/s)'); ylabel('Imaginary part (1/s)');
    elseif n==3
      ids=1:4; hh=[]; labels={};
      for j=ids, hh(end+1)=plot(times{j},responses{j},'linewidth',1.8); labels{end+1}=sprintf('K=%.5g',K(j)); end
      plot([0 55],[.98 .98],'k:'); plot([0 55],[1.02 1.02],'k:');
      xlim([0 55]); xlabel('Time (s)'); ylabel('Unit-step response'); title('C: Stable gain comparisons; dotted lines = 2% band'); legend(hh,labels,'location','northeast');
    elseif n==4
      for j=[5 6], plot(times{j},responses{j},'linewidth',1.6); end
      xlim([0 80]); xlabel('Time (s)'); ylabel('Unit-step response'); title('C: Boundary oscillation and unstable growth'); legend('K=48: sustained oscillation','K=50: growing oscillation','location','northwest');
    else
      t=times{2}; y2=step(tf(16/9,[1 4/3 16/9]),t);
      h1=plot(t,responses{2},'b','linewidth',2); h2=plot(t,y2,'r--','linewidth',1.8);
      plot([0 12],[.98 .98],'k:'); plot([0 12],[1.02 1.02],'k:'); xlim([0 12]);
      xlabel('Time (s)'); ylabel('Unit-step response'); title(sprintf('B: Selected K=%.6f, dominant zeta=0.5',Ksel));
      legend([h1 h2],{'Full third-order model','Dominant second-order prediction'},'location','northeast');
      text(6,.65,sprintf('Actual Mp=%.3f%%, ts=%.3f s',metrics(2,1),metrics(2,4)));
      text(6,.55,sprintf('Predicted Mp=%.3f%%, ts approx %.3f s',predicted(2,1),predicted(2,2)));
    end
    grid on; box on; set(gca,'fontsize',12);
    print(f,fullfile(out,[stems{n} '.png']),'-dpng','-r150'); print(f,fullfile(out,[stems{n} '.svg']),'-dsvg'); close(f);
  end
end
