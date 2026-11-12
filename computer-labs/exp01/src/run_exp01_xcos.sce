// Execute with Scilab 2026.1 in console mode with JVM support.
// Local launcher must redirect SCIHOME and TEMP into the project .local/.
src=get_absolute_file_path("run_exp01_xcos.sce");
root=fullpath(src+"../../");
out=root+".local/exp01/results/";
if ~isdir(out) then mkdir(out); end
loadXcosLibs(); loadScicos();
exec(src+"build_exp01.sci",-1);
diary(out+"xcos_run.txt");
mprintf("Scilab %s, native Xcos fixed-step RK4, zero initial state\n",getversion());
for K=[1,4,10]
    for h=[0.02,0.01,0.005]
        scs_m=build_exp01(K,h);
        filename=out+msprintf("feedback_K%d_h%g.zcos",K,h);
        xcosDiagramToScilab(filename,scs_m);
        // Reimport the saved editable model before running it.
        scs_m=xcosDiagramToScilab(filename);
        clear exp01_scope;
        Info=scicos_simulate(scs_m,list(),struct(),"nw");
        times=matrix(exp01_scope.time,-1,1);
        values=matrix(exp01_scope.values,-1,1);
        wd=sqrt(4+K);
        exact=K/(5+K)*(1-exp(-times).*(cos(wd*times)+sin(wd*times)/wd));
        err=max(abs(values-exact));
        if err>1d-5 then error("Xcos/analytic discrepancy exceeded tolerance"); end
        csvWrite([times,values],out+msprintf("xcos_K%d_h%g.csv",K,h),",",".","%.17g",["time_s,displacement_m"]);
        mprintf("K=%g h=%g: %d samples, max analytic error %.12g\n",K,h,length(times),err);
        if K==4 & h==0.01 then
            xcosDiagramExport(scs_m,out+"05_xcos_model.png");
            xcosDiagramExport(scs_m,out+"05_xcos_model.svg");
        end
    end
end
diary(0);
exit;
