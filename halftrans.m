function [t]= halftrans(i,j)

global A visc dx perm
if perm(i) > 0
 t=((A/(visc*dx))*2/(1/perm(j)+1/perm(i)));
else 
    t=0;
end
end
