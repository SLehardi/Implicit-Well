clear all
clc

global Nx Ny N dx dy x y h A L W V perm phi visc cf dt Pinit J req rw s
myinputs;

trans = (perm*A/(visc*dx));       % Calculate transmissibility

% Set up T and B arrays and vector Q
T=zeros(N,N);B=zeros(N,N); Q=zeros(N,1); 

for k=1:N % iterates for all the rows in the matrix
   
   Tsum=0;
   if rem(k,Nx)~=0 % if k (row number) is NOT divisible by m, then set t(k,k+1) = -3
   T(k,k+1)=-halftrans(k,k+1);
   Tsum= Tsum-T(k,k+1);
   end
   if rem(k,Nx)~=1 % if the remainder is NOT equal to 1, this means that we are NOT in the left bound, then the element is -3
   T(k,k-1)=-halftrans(k,k-1);
   Tsum= Tsum-T(k,k-1);
   end
   if (k/Nx)>1 % if the division of k by m is higher than 1, we are NOT at the bottom bound
   T(k,k-Nx)=-halftrans(k,k-Nx);
   Tsum= Tsum-T(k,k-Nx);
   end
   if ceil(k/Nx)<Ny % if k/m is NOT rounded up to n, we are NOT at the upper bound
   T(k,k+Nx)=-halftrans(k,k+Nx);
   Tsum= Tsum-T(k,k+Nx);
   end
   T(k,k)= Tsum;
   B(k,k)=V*phi(k)*cf; 
end

 [Jmatrix,Q] = Wells;


time =0; c=1;
AA=(T*6.33*10^-3+Jmatrix+B/dt);
AA=sparse(AA);

Pold=Pinit*ones(N,1);
while time < 10
    
   
    BB=((1/dt)*B*Pold+Q);
    Pnew= AA\BB;
    Pold=Pnew;
    P_plot(:,c)=Pnew;
    time=time +dt;
    Time(c)=time;
    qw(c)=J*(Pnew(2604)-500);
    c=c+1;
    
end    

x1=find(perm < 0.01);
P_plot(x1,:)=NaN;
[X,Y]= meshgrid(x,y);
P_early= P_plot(:,1);
P_mid= P_plot(:,10);
P_late=P_plot(:,100);
P_early = reshape(P_early, [80 60])' ;
P_mid = reshape(P_mid, [80 60])' ;
P_late = reshape(P_late, [80 60])' ;

figure %open a new figure
contourf(X,Y,P_early,100) %creates filled contour map
title('Pressure at Early Time') % Puts a title on the figure
xlabel('X') %Adds X axis label
ylabel('Y') %Adds Y axis label
colormap('jet') %Changes the color scheme to warm colors
colorbar %Adds colorbar to plot

figure %open a new figure
contourf(X,Y,P_mid,100) %creates filled contour map
title('Pressure at Intermediate Time') % Puts a title on the figure
xlabel('X') %Adds X axis label
ylabel('Y') %Adds Y axis label
colormap('jet') %Changes the color scheme to warm colors
colorbar %Adds colorbar to plot

figure %open a new figure
contourf(X,Y,P_late,100) %creates filled contour map
title('Pressure at Late time') % Puts a title on the figure
xlabel('X') %Adds X axis label
ylabel('Y') %Adds Y axis label
colormap('jet') %Changes the color scheme to warm colors
colorbar %Adds colorbar to plot

%Bonus 1
figure
plot(Time, qw);
title('Flowrate vs. Time') % Puts a title on the figure
xlabel('Time') %Adds X axis label
ylabel('Flowrate ft3/day') %Adds Y axis label

%Bonus 2
figure
perm_plot=perm;
perm_plot(x1,:)=NaN;
perm_plot = reshape(perm_plot, [80 60])' ;
contourf(X,Y, perm_plot)
title('Permeability') % Puts a title on the figure
xlabel('X') %Adds X axis label
ylabel('Y') %Adds Y axis label
colorbar %Adds colorbar to plot

