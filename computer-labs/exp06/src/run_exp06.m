function run_exp06(mode = 'all')
  % Experiment 06. Frequency in rad/s, magnitude in dB, phase in degrees.
  % GNU Octave 11.3 + control 4.2.1. Outputs are private and reproducible.
  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  out = fullfile(root, '.local', 'exp06', 'results');
  if (!exist(out, 'dir')), error('Prepare private low-integrity results directory first.'); endif
  pkg load control;
  if (any(strcmp(mode, {'all','compute'}))), compute_ab(out); endif
  if (any(strcmp(mode, {'all','fit'}))), compute_c(out); endif
  if (!any(strcmp(mode, {'all','compute','fit'}))), error('mode: all, compute, fit'); endif
endfunction

function [m,p] = truth(w)
  h = 10 ./ ((1i*w) .* (1i*w+2) .* (1i*w+5));
  m = 20*log10(abs(h));
  p = -90 - atan(w/2)*180/pi - atan(w/5)*180/pi;
endfunction

function [m,p,a2,a5] = asym(w)
  m = -20*log10(w)-20*max(0,log10(w/2))-20*max(0,log10(w/5));
  a2 = -45*min(2,max(0,log10(w/2)+1));
  a5 = -45*min(2,max(0,log10(w/5)+1));
  p = -90+a2+a5;
endfunction

function write_csv(out,name,header,x)
  f=fopen(fullfile(out,name),'w'); if(f<0),error('Cannot write CSV');endif
  fprintf(f,'%s\n',header); fclose(f);
  dlmwrite(fullfile(out,name),x,'-append','delimiter',',','precision','%.16g');
endfunction

function compute_ab(out)
  w=unique([logspace(-1,2,6001),.2,.5,1,2,sqrt(10),5,10,20,50])';
  [m,p]=truth(w); [a,b,a2,a5]=asym(w);
  write_csv(out,'a_factors.csv',strjoin({'w_rad_s','integrator_mag_db','pole2_mag_db','pole5_mag_db','exact_mag_db','asym_mag_db','integrator_phase_deg','pole2_phase_deg','pole5_phase_deg','exact_phase_deg','asym_pole2_phase_deg','asym_pole5_phase_deg','asym_phase_deg'},','),[w,-20*log10(w),-10*log10(1+(w/2).^2),-10*log10(1+(w/5).^2),m,a,-90*ones(size(w)),-atan(w/2)*180/pi,-atan(w/5)*180/pi,p,a2,a5,b]);
  header='w_rad_s,exact_mag_db,asym_mag_db,mag_error_db,exact_phase_deg,asym_phase_deg,phase_error_deg';
  write_csv(out,'b_dense.csv',header,[w,m,a,a-m,p,b,b-p]);
  wp=[.2;1;2;5;10;50]; [mp,pp]=truth(wp);[ap,bp]=asym(wp);
  write_csv(out,'b_points.csv',header,[wp,mp,ap,ap-mp,pp,bp,bp-pp]);
  % Exhaustive smooth-interval stationary points + all phase-fold endpoints.
  % d(error)/d ln(w) = -45*n/ln(10) + (180/pi)*(2w/(4+w^2)+5w/(25+w^2)).
  edges=[.1,.2,.5,2,5,20,50,100]; candidates=edges;
  for j=1:numel(edges)-1
    mid=sqrt(edges(j)*edges(j+1)); active=(mid>.2 && mid<20)+(mid>.5 && mid<50);
    d=45*active/log(10);
    if(active>0)
      r=roots([d,-7*180/pi,29*d,-70*180/pi,100*d]);
      r=real(r(abs(imag(r))<1e-9 & real(r)>edges(j) & real(r)<edges(j+1)));
      candidates=[candidates,r'];
    endif
  endfor
  pe=arrayfun(@phase_error,candidates); [pem,idx]=max(pe);
  [mm,pm]=truth([2;5]); [am,bm]=asym([2;5]);
  % Preserve all tied maxima rather than selecting floating-point tail noise.
  tied=edges(arrayfun(@phase_error,edges)>=pem-1e-10);
  maxima=[1,2,am(1)-mm(1),1;1,5,am(2)-mm(2),1];
  for j=1:numel(tied), maxima(end+1,:)=[2,tied(j),phase_error(tied(j)),2]; endfor
  assert(numel(tied)==2 && all(abs(tied-[.2,50])<1e-10));
  assert(abs(max(abs(b-p))-pem)<1e-10);
  write_csv(out,'b_maxima.csv','quantity_code,w_rad_s,max_abs_error,unit_code',maxima);
  % Scalar control checks avoid wide-vector polynomial fitting in freqresp.
  s=tf('s'); G=10/(s*(s+2)*(s+5)); rel=0;
  for k=1:numel(wp)
    hc=freqresp(G,wp(k)); hd=10/((1i*wp(k))*(1i*wp(k)+2)*(1i*wp(k)+5));
    rel=max(rel,abs(hc-hd)/abs(hd));
  endfor
  assert(rel<1e-11); assert(abs((am(1)-mm(1))-(am(2)-mm(2)))<1e-12);
  assert(abs(max(a-m)-(am(1)-mm(1)))<1e-11);
  f=fopen(fullfile(out,'ab_log.txt'),'w');
  fprintf(f,'Octave %s; control scalar freqresp relative error %.16g\n',version,rel);
  fprintf(f,'Low-frequency 0.2 rad/s baseline %.12f dB\n',-20*log10(.2));
  fprintf(f,'Amplitude maxima at BOTH 2 and 5 rad/s: %.12f dB\n',am(1)-mm(1));
  fprintf(f,'Phase global maximum %.12f deg at BOTH 0.2 and 50 rad/s\n',pem);
  fclose(f); disp(fileread(fullfile(out,'ab_log.txt')));
endfunction

function e=phase_error(w)
  [m,p]=truth(w);[a,b]=asym(w);e=abs(b-p);
endfunction

function compute_c(out)
  w=[.1;.16;.25;.4;.63;1;1.6;2.5;4;6.3;10;16;25;40;63;100];
  [m,p]=truth(w); mr=round(m*10)/10; pr=round(p);
  subset=[1,3,5,7,9,11,13,16];
  h='w_rad_s,mag_db,phase_deg';
  write_csv(out,'c_input_full16.csv',h,[w,m,p]);
  write_csv(out,'c_input_rounded16.csv',h,[w,mr,pr]);
  write_csv(out,'c_input_rounded8.csv',h,[w(subset),mr(subset),pr(subset)]);
  % The estimator receives DATA ONLY, never the true model coefficients.
  [fit1,coarse,diag1,starts1]=identify_exp06(w,m,p);
  [fit2,c2,diag2,starts2]=identify_exp06(w,mr,pr);
  [fit3,c3,diag3,starts3]=identify_exp06(w(subset),mr(subset),pr(subset));
  models=[coarse;fit1;fit2;fit3];
  % Half-grid shift and remove input coincidences: validation points held out.
  v=10.^(-1+((0:7999)'+.5)*3/8000);
  for j=1:numel(w), v(abs(v-w(j))<1e-12)=[]; endfor
  [mt,pt]=truth(v); vals=[v,mt,pt]; errors=zeros(4,4);
  for j=1:4
    [me,ph]=estimate_response(models(j,:),v);vals=[vals,me,ph];
    [em,im]=max(abs(me-mt));[ep,ip]=max(abs(ph-pt)); errors(j,:)=[em,v(im),ep,v(ip)];
  endfor
  write_csv(out,'c_validation.csv','w_rad_s,true_mag_db,true_phase_deg,coarse_mag_db,coarse_phase_deg,full16_mag_db,full16_phase_deg,rounded16_mag_db,rounded16_phase_deg,rounded8_mag_db,rounded8_phase_deg',vals);
  % Model = K/[s*(1+s/p1)*(1+s/p2)] = numerator/[s*(s+p1)*(s+p2)].
  write_csv(out,'c_models.csv','model_code,K,p1_rad_s,p2_rad_s,numerator,den_s2,den_s1,max_mag_error_db,max_mag_error_w,max_phase_error_deg,max_phase_error_w',[ (1:4)',models,prod(models,2),sum(models(:,2:3),2),prod(models(:,2:3),2),errors]);
  write_csv(out,'c_identification.csv','case_code,n_points,low_slope_db_dec,high_slope_db_dec,origin_poles,finite_poles,coarse_K,coarse_p1,coarse_p2,objective',[1,diag1;2,diag2;3,diag3]);
  write_csv(out,'c_multistart.csv','case_code,start_K,start_p1,start_p2,fit_K,fit_p1,fit_p2,objective,exitflag',[ones(rows(starts1),1),starts1;2*ones(rows(starts2),1),starts2;3*ones(rows(starts3),1),starts3]);
  assert(errors(2,1)<1e-7 && errors(2,3)<1e-6);
  assert(all(isfinite(vals(:))) && all(models(:)>0));
  f=fopen(fullfile(out,'c_log.txt'),'w');
  fprintf(f,'Data: generated from the course model, not physical measurements.\n');
  fprintf(f,'Model codes: 1 coarse full16; 2 fitted full16; 3 fitted rounded16; 4 fitted rounded8.\n');
  fprintf(f,'Rounded precision 0.1 dB and 1 degree; rounded8 subset indices: 1 3 5 7 9 11 13 16.\n');
  fprintf(f,'Independent validation: %d held-out frequencies, 0.1 to 100 rad/s.\n',numel(v));
  for j=1:4
    fprintf(f,'model %d: K=%.12g p1=%.12g p2=%.12g maxerr=%.12g dB, %.12g deg\n',j,models(j,:),errors(j,1),errors(j,3));
  endfor
  fclose(f);disp(fileread(fullfile(out,'c_log.txt')));
endfunction

function [m,p]=estimate_response(q,w)
  m=20*log10(q(1))-20*log10(w)-10*log10(1+(w/q(2)).^2)-10*log10(1+(w/q(3)).^2);
  p=-90-atan(w/q(2))*180/pi-atan(w/q(3))*180/pi;
endfunction
