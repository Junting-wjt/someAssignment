%--------------------------------------------
%    Rectangular Membrane Mode Plotter 
%           Dr Thomas McKenzie
%           Reid School of Music 
%          University of Edinburgh               
%--------------------------------------------

clear all
close all

%-- Global Parameters

Lx      = 1;
Ly      = 0.75;
N       = 4; 

%-- derived paramters
res     = 0.02;
dx      = res*Lx;
dy      = res*Ly;

x       = dx : dx : Lx - dx;
y       = dy : dy : Ly - dy;

omMat      = zeros(N^2,3);
ind        = 0;

for m = 1 : N
    for n = 1 : N
        ind             = ind + 1;
        omMat(ind,1)    = sqrt(m^2/Lx^2 + n^2/Ly^2); 
        omMat(ind,2)    = m;
        omMat(ind,3)    = n;
    end
end

%% -- plot 

[X,Y]   = meshgrid(x,y);
v       = [0,0];

figure; 
for q = 1 : 12
    m       = omMat(q,2);
    n       = omMat(q,3);
    w       = sin(m*X*pi/Lx) .* sin(n*Y*pi/Ly);
    subplot(3,4,q);
    mesh(X,Y,w,'FaceAlpha',0.5); hold on; %,'FaceColor','flat' %'FaceLighting','flat','EdgeLighting','flat',
        title(['m=',num2str(m),', n=',num2str(n)])
    view([90,20,85]);
    contour(X, Y, w,v,'color','k','linewidth',2)
end