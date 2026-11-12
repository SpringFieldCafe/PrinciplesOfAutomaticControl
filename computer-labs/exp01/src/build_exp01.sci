// Construct an editable native Xcos feedback model from library blocks.
// Units: m=1 kg, b=2 N*s/m, k=5 N/m; K N/m; time in s.
// Each block is an independently editable Xcos component.
function diagram = build_exp01(K, h)
    diagram=scicos_diagram();
    diagram.props.title="Experiment 01 - mass spring damper feedback";
    diagram.props.tf=10;
    diagram.props.context=["K="+string(K)+";";"h="+string(h)+";"];
    diagram.props.tol(1)=1d-9; diagram.props.tol(2)=1d-9;
    diagram.props.tol(6)=6; // Native fixed-step Runge-Kutta 4 solver.
    diagram.props.tol(7)=h;

    step=STEP_FUNCTION("define");
    step.model.rpar.objs(1).model.firing=0;
    step.model.rpar.objs(1).model.rpar=[1;1];
    step.model.rpar.objs(1).graphics.exprs=["0";"0";"1"];
    step.graphics.orig=[30,180]; step.graphics.sz=[50,40];
    step.graphics.id="Unit reference (m)";
    diagram.objs(1)=step;

    sumblock=SUMMATION("define");
    sumblock.graphics.orig=[150,175]; sumblock.graphics.sz=[40,50];
    sumblock.graphics.exprs=["1";"[1;-1]";"0"];
    diagram.objs(2)=sumblock;

    gain=GAINBLK_f("define"); gain.model.rpar=K;
    gain.graphics.exprs="K"; gain.graphics.orig=[250,180]; gain.graphics.sz=[60,40];
    diagram.objs(3)=gain;

    plant=CLR("define");
    A=[0,1;-5,-2]; B=[0;1]; C=[1,0]; D=0;
    plant.model.state=[0;0]; plant.model.rpar=[A(:);B(:);C(:);D];
    plant.graphics.exprs=["1";"s^2+2*s+5"];
    plant.graphics.orig=[370,175]; plant.graphics.sz=[120,50];
    plant.graphics.id="Displacement plant";
    diagram.objs(4)=plant;

    split=SPLIT_f("define"); split.graphics.orig=[530,200];
    diagram.objs(5)=split;
    sink=TOWS_c("define"); name="exp01_scope";
    sink.model.ipar=[20005;length(name);ascii(name)'];
    sink.graphics.exprs=["20005";name;"0"];
    sink.graphics.orig=[620,240]; sink.graphics.sz=[80,40];
    diagram.objs(6)=sink;
    scope=CSCOPE("define"); scope.model.rpar=[0;-0.1;1.1;10];
    scope.graphics.exprs(5)="-0.1"; scope.graphics.exprs(6)="1.1";
    scope.graphics.exprs(7)="10";
    scope.graphics.orig=[620,140]; scope.graphics.sz=[80,50];
    diagram.objs(7)=scope;
    clock=CLOCK_c("define");
    clock.model.rpar.objs(2).model.rpar=[h;0];
    clock.model.rpar.objs(2).model.firing=0;
    clock.model.rpar.objs(2).graphics.exprs=["h";"0"];
    clock.graphics.orig=[600,350]; clock.graphics.sz=[50,40];
    diagram.objs(8)=clock;
    splitclock=CLKSPLIT_f("define"); splitclock.graphics.orig=[660,320];
    diagram.objs(9)=splitclock;

    diagram=wire(diagram,1,1,2,1,[80;150],[200;208],%f);
    diagram=wire(diagram,2,1,3,1,[190;250],[200;200],%f);
    diagram=wire(diagram,3,1,4,1,[310;370],[200;200],%f);
    diagram=wire(diagram,4,1,5,1,[490;530],[200;200],%f);
    diagram=wire(diagram,5,1,2,2,[530;530;120;120;150],[200;80;80;192;192],%f);
    diagram=wire(diagram,5,2,6,1,[530;570;570;620],[200;200;260;260],%f);
    diagram=wire(diagram,5,3,7,1,[530;570;570;620],[200;200;165;165],%f);
    diagram=wire(diagram,8,1,9,1,[625;625;660],[350;320;320],%t);
    diagram=wire(diagram,9,1,6,1,[660;660],[320;280],%t);
    diagram=wire(diagram,9,2,7,1,[660;740;740;660;660],[320;320;210;210;190],%t);
endfunction

function diagram=wire(diagram,src,sp,dst,dp,xx,yy,event)
    n=length(diagram.objs)+1;
    link=scicos_link(); link.from=[src,sp,0]; link.to=[dst,dp,1];
    link.xx=xx; link.yy=yy;
    if event then
        link.ct=[5,-1];
        diagram.objs(src).graphics.peout(sp)=n;
        diagram.objs(dst).graphics.pein(dp)=n;
    else
        diagram.objs(src).graphics.pout(sp)=n;
        diagram.objs(dst).graphics.pin(dp)=n;
    end
    diagram.objs(n)=link;
endfunction
