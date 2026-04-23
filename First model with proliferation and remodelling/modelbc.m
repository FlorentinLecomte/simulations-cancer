function [pl,ql,pr,qr] = modelbc(xl,ul,xr,ur,t) 
pl = [0; 0; 0];
ql = [-1;1; -1];
pr = pl;
qr = ql;
end