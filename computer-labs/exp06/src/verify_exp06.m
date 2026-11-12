function verify_exp06()
  root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out=fullfile(root,'.local','exp06','results');
  c=dlmread(fullfile(out,'c_models.csv'),',',1,0);
  v=dlmread(fullfile(out,'c_validation.csv'),',',1,0); w=v(:,1);err=0;
  for j=1:4
    h=c(j,5)./polyval([1,c(j,6),c(j,7),0],1i*w);
    mm=20*log10(abs(h));pp=unwrap(angle(h))*180/pi;
    err=max(err,max(abs(mm-v(:,2*j+2))));
    err=max(err,max(abs(pp-v(:,2*j+3))));
  endfor
  assert(err<1e-10);
  % Different synthetic parameters test data-only behavior of the estimator.
  ww=logspace(-2,3,18)';
  mm=20*log10(1.4)-20*log10(ww)-10*log10(1+(ww/1.3).^2)-10*log10(1+(ww/7.2).^2);
  pp=-90-atan(ww/1.3)*180/pi-atan(ww/7.2)*180/pi;
  q=identify_exp06(ww,mm,pp);e=max(abs(q-[1.4,1.3,7.2]));assert(e<1e-7);
  f=fopen(fullfile(out,'independent_checks.txt'),'w');
  fprintf(f,'Expanded polynomial vs factor-form CSV max absolute difference %.16g (dB or degrees).\n',err);
  fprintf(f,'Data-only estimator on separate synthetic K=1.4,p1=1.3,p2=7.2: fitted %.14g %.14g %.14g; max parameter error %.16g.\n',q,e);
  fclose(f);disp(fileread(fullfile(out,'independent_checks.txt')));
endfunction
