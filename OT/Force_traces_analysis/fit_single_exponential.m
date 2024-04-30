function [f] = fit_single_exponential(ff,ee)
    % double exponential
    %p0=[ee(1) 0.1 5 1]; lb=[0 0 0 0]; ub=[10 10 25 10];
    %ft = fittype('e0.*exp(x./f0)+e1.*exp(x./f1)');
    %single exponential
    p0 = [ee(1) 5]; lb = [0 0]; ub = [10 25];
    ft = fittype('e0.*exp(x./f0)');
    
    options = fitoptions(ft);
    options.startPoint = p0;
    options.Lower = lb;
    options.Upper = ub;
    jkeep = find(isnan(ee)==0 & ee ~= 0);
    f = fit(ff(jkeep)',ee(jkeep)',ft,options);
end