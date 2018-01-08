function myinputs
global Nx Ny N dx dy x y h A L W V perm phi visc cf dt Pinit req rw s
% Define reservoir properties and boundary conditions (Need to be entered)
L=4000;             % Reservoir length (ft)
W= 3000;
h=30;
Nx=80;
Ny=60;
N=Nx*Ny;               % Number of grids
dx = L/Nx;                    % Grid size
dy= W/Ny;
x=linspace(dx/2,L-dx/2,Nx);   % Postion of grid centers
y=linspace(dy/2,W-dy/2,Ny);
req= .20*dx;
rw= 4/12;
s=0;
A=dx*h;            % Cross Sectional Area (ft2)
V= dx*dy*h;
%perm=30;               % Permeability (md)
fid=fopen('permeability.txt');
perm=fscanf(fid,'%f');
%phi=.25;             % Porosity 
fid=fopen('porosity.txt');
phi=fscanf(fid,'%f');
visc=1;             % Viscosity (cp)
cf=10^-6;           % Compressibility (psi^-1)
dt=0.1;            % Timestep (days)
Pinit=1000;   % Initial Reservoir Pressure

