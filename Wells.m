function [Jmatrix,Q] = Wells

 global J Nx Ny dx dy perm visc req rw s h
% well inputs containing x,y, type of well, and pressure
 Well= [ 855 1102 1 1300;...
        3215 800 1 1300;...
        1004 1992 1 1300;...
        3100 2450 1 1300;...
        2155 1650 2 500;]';

Nw= length(Well);
Pro= Nx*Ny;
Jmatrix= zeros(Pro,Pro); Q=zeros(Pro,1);J= zeros(Pro,1);

for N = 1:Nw
    % finding ratios of x and y to determine well position in l coordinate
    % system
    i= ceil(Well(1,N)/dx);
    j= ceil(Well(2,N)/dy);
    l=(j-1)*Nx+i;
    J= 2*6.33*10^-3*pi*perm(l)*h/(visc*log(req/rw)+s);
    Type=Well(3,N);
    % Determines what type of well and builds Q and J matrix
    if Type ==1
        Q(l,1)= Well(4,N);
    
    elseif Type==0
         Q(l,1)= -Well(4,N);
    
    else
        Q(l,1)= J*Well(4,N);
        Jmatrix(l,l)= J;
    end
            
            
end
        
        


