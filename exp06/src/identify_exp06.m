function [best,coarse,diagnostic,starts] = identify_exp06(w,mag_db,phase_deg)
  % Data-only estimator; assumes positive gain, no zeros, real stable poles.
  % Infer integrator and finite pole count from endpoint magnitude slopes.
  w=w(:);mag_db=mag_db(:);phase_deg=phase_deg(:);
  assert(numel(w)>=8 && all(diff(w)>0) && all(w>0));
  assert(all(isfinite([w;mag_db;phase_deg])));
  x=log10(w); n=numel(w); k=min(3,floor(n/3));
  lo=polyfit(x(1:k),mag_db(1:k),1); hi=polyfit(x(end-k+1:end),mag_db(end-k+1:end),1);
  n0=round(-lo(1)/20);nf=round(-hi(1)/20)-n0;
  assert(n0==1 && nf==2,'Data endpoint slopes do not support one integrator + two finite poles.');
  slopes=diff(mag_db)./diff(x); xm=(x(1:end-1)+x(2:end))/2;
  % Locate -30/-50 dB/dec transitions by interpolation in local slope.
  p1=10^crossing(xm,slopes,-30);p2=10^crossing(xm,slopes,-50);
  K=10^(mean(mag_db(1:k)+20*x(1:k))/20);
  coarse=[K,p1,p2];
  % All starts derive from coarse data estimates, with generic scale offsets.
  scales=[.6,1,1.7];starts=[];best=[];bestcost=Inf;
  opts=optimset('Display','off','MaxIter',1500,'MaxFunEvals',4000,'TolX',1e-11,'TolFun',1e-16);
  for a=scales
    for b=scales
      initial=[K,p1*a,p2*b]; initial(2:3)=sort(initial(2:3));
      initial(3)=max(initial(3),initial(2)*1.05);
      theta=log([initial(1),initial(2),initial(3)-initial(2)]);
      [th,cost,flag]=fminsearch(@(t)objective(t,w,mag_db,phase_deg),theta,opts);
      q=decode(th); starts(end+1,:)=[initial,q,cost,flag];
      if(cost<bestcost),bestcost=cost;best=q;endif
    endfor
  endfor
  assert(isfinite(bestcost) && any(starts(:,8)>0),'No converged fit');
  diagnostic=[n,lo(1),hi(1),n0,nf,coarse,bestcost];
endfunction

function z=crossing(x,s,target)
  j=find(s(1:end-1)>=target & s(2:end)<=target,1);
  if(isempty(j)),[v,j]=min(abs(s-target));z=x(j);
  else,z=x(j)+(target-s(j))/(s(j+1)-s(j))*(x(j+1)-x(j));endif
endfunction

function q=decode(t)
  e=exp(t);q=[e(1),e(2),e(2)+e(3)];
endfunction

function cost=objective(t,w,m,p)
  if(any(abs(t)>50)),cost=1e100;return;endif
  q=decode(t);
  mm=20*log10(q(1))-20*log10(w)-10*log10(1+(w/q(2)).^2)-10*log10(1+(w/q(3)).^2);
  pp=-90-atan(w/q(2))*180/pi-atan(w/q(3))*180/pi;
  % Explicit weights: one dB and ten degrees contribute equal squared cost.
  cost=mean((mm-m).^2+((pp-p)/10).^2);
endfunction
