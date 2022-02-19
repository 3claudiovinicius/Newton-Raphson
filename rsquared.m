function [R2, theta] = rsquared(Hist, par)
Vhist = Hist.Values;
theta = Hist.BinCounts(1)/length(Hist.Data);
a2 = wblcdf(Hist.BinEdges(find(Hist.BinEdges>0)), par(1), par(2));
a1 = wblcdf(Hist.BinEdges(find(Hist.BinEdges<15)), par(1), par(2));

Vest = (a2-a1)*(1-theta)/(Hist.BinWidth); %Valores estimados pelos parâmetros alpha e beta
Vest(1) = (theta/Hist.BinWidth)+Vest(1);
r = corrcoef(Vhist, Vest);
R2 = r(2)^2;



end