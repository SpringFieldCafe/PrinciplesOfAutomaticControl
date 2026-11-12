function exp02_diagram()
% Editable signal-flow diagram; executable equations are in run_exp02.m.
  root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out=fullfile(root,'.local','exp02','results');
  if ~exist(out,'dir'),mkdir(out);end
  f=figure('visible','off','position',[50 50 1200 550]); axis([0 12 0 5]); axis off; hold on;
  box(3,2.7,1.5,0.8,'C(s) = K'); box(7,2.7,2.5,0.8,'G(s) = 1/[s(s+2)]');
  text(0.1,3.1,'r','fontsize',14); arrow(0.5,3.1,1.6,3.1); circle(2,3.1); text(1.85,3.1,'+','fontsize',14);
  arrow(2.25,3.1,3,3.1); text(2.45,3.4,'e'); arrow(4.5,3.1,5.55,3.1); text(4.9,3.4,'u');
  circle(5.8,3.1);text(5.67,3.1,'+','fontsize',14);arrow(6.05,3.1,7,3.1);
  text(5.7,4.35,'d'); arrow(5.8,4.1,5.8,3.35); text(6,3.7,'+');
  arrow(9.5,3.1,11.2,3.1);text(11.4,3.1,'y','fontsize',14);
  plot([10.3 10.3],[3.1 1.2],'k-'); arrow(10.3,1.2,8.25,1.2);circle(8,1.2);text(7.87,1.2,'+','fontsize',14);
  text(8,0.2,'n');arrow(8,0.5,8,0.95); text(8.2,0.7,'+');
  plot([7.75 2],[1.2 1.2],'k-');arrow(2,1.2,2,2.85);text(1.6,2.5,'-','fontsize',18);
  text(4.2,0.9,'measured output y + n','fontsize',12);
  text(0.3,4.7,'Negative unity feedback: plant-input load d and measurement noise n','fontsize',15);
  print(f,fullfile(out,'signal_flow.png'),'-dpng','-r150');print(f,fullfile(out,'signal_flow.svg'),'-dsvg');close(f);
end
function box(x,y,w,h,label)
  rectangle('position',[x y w h],'linewidth',1.6);text(x+w/2,y+h/2,label,'horizontalalignment','center','fontsize',12);
end
function circle(x,y)
  rectangle('position',[x-.25 y-.25 .5 .5],'curvature',[1 1],'linewidth',1.5);
end
function arrow(x1,y1,x2,y2)
  quiver(x1,y1,x2-x1,y2-y1,0,'k','linewidth',1.5,'maxheadsize',.15);
end
