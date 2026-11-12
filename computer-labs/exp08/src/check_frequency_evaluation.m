function check_frequency_evaluation()
% Isolate vector freqresp conditioning from direct Horner evaluation.
  pkg load control;
  root=fileparts(fileparts(fileparts(mfilename('fullpath'))));out=fullfile(root,'.local','exp08','results');
  s=tf('s');L=48*(s+2)/(s+12)*(s+.05)/(s+1/30)/(s*(s+2)*(s+5));
  S=feedback(tf(1),L);T=feedback(L,1);w=logspace(-4,4,5001)';
  [ns,ds]=tfdata(S,'v');[nt,dt]=tfdata(T,'v');
  sv=squeeze(freqresp(S,w));tv=squeeze(freqresp(T,w));
  sh=polyval(ns,1i*w)./polyval(ds,1i*w);th=polyval(nt,1i*w)./polyval(dt,1i*w);
  sampled=unique(round(linspace(1,numel(w),51)));rows=[];
  for ix=sampled
    scalar=squeeze(freqresp(T,w(ix)));direct=polyval(nt,1i*w(ix))/polyval(dt,1i*w(ix));
    rows=[rows;w(ix),abs(tv(ix)-th(ix)),abs(scalar-th(ix)),abs(direct-th(ix))];
  end
  assert(max(rows(:,3))<1e-10);assert(max(rows(:,4))<1e-10);
  fid=fopen(fullfile(out,'frequency_diagnostic.txt'),'w');
  fprintf(fid,'control 4.2.1 @tf/__freqresp__ uses __polyval__ + polyfit when vector length exceeds polynomial degree. Source inspected read-only.\n');
  fprintf(fid,'wide vector freqresp max |S+T-1| = %.17g\n',max(abs(sv+tv-1)));
  fprintf(fid,'direct vector polyval max |S+T-1| = %.17g\n',max(abs(sh+th-1)));
  fprintf(fid,'51 scalar freqresp vs direct vector Horner max error = %.17g\n',max(rows(:,3)));
  fprintf(fid,'51 scalar vs vector polyval max error = %.17g\n',max(rows(:,4)));
  fprintf(fid,'Fault isolated to control vector frequency evaluation / polynomial refit route, not direct vector polyval. Installed software unchanged.\n');fclose(fid);
  path=fullfile(out,'frequency_diagnostic.csv');fid=fopen(path,'w');fprintf(fid,'omega,wide_freqresp_error,scalar_freqresp_error,scalar_polyval_error\n');fclose(fid);dlmwrite(path,rows,'-append','delimiter',',','precision','%.17g');
  disp(fileread(fullfile(out,'frequency_diagnostic.txt')));
  load(fullfile(out,'models.mat'),'Cs','G');checks=[];
  for j=1:numel(Cs)
    L=Cs{j}*G;[gm,pm,wp,wg]=margin(L);[n,d]=tfdata(L,'v');
    hg=polyval(n,1i*wg)/polyval(d,1i*wg);hp=polyval(n,1i*wp)/polyval(d,1i*wp);
    errors=[abs(abs(hg)-1),abs((180+angle(hg)*180/pi)-pm),abs(gm*abs(hp)-1),abs(imag(hp))];
    assert(max(errors)<1e-9);checks=[checks;j-1,errors];
  end
  path=fullfile(out,'margin_crosscheck.csv');fid=fopen(path,'w');fprintf(fid,'id,gain_cross_error,PM_error_deg,GM_error,phase_cross_imag_error\n');fclose(fid);dlmwrite(path,checks,'-append','delimiter',',','precision','%.17g');
end
