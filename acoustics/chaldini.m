%-----------------------------------------------------
%            Chladni Patterns         
%           Dr Thomas McKenzie
%           Reid School of Music 
%          University of Edinburgh   
%    Plotting degenerate mode shapes of a 2D membrane 
%-----------------------------------------------------

clc
clear all
close all

% input params
L = 1;              % membrane size (square)
Nx = 100;           % num grid points in x,y
m = 3; n = 2;       % modal indices
bc = 0;             % choice of boundary condition -> 0: free; 1: fixed

% weighting of modes (for chladni patterns)
%theta = 7*pi/4;       % a = cos(theta); b = sin(theta)

% time-varying theta - uncomment for animation
Nt = 100;                    % num frames in animation of one cycle
theta = 2*pi*(0:Nt-1)/Nt;    % discrete angles

% Mode shapes
if bc == 0
    % free boundary conditions - centered around x=y=0 
    x = linspace(-L/2, L/2, Nx);             % discretize x axis domain
    [X,Y] = meshgrid(x,x);                  % create 2D grid of points
    Wmn = cos(m*pi*X/L) .* cos(n*pi*Y/L);    % mode shapes
    Wnm = cos(n*pi*X/L) .* cos(m*pi*Y/L);
else
    % fixed boundary conditions - fixed at x=y=0 and x=y=L
    x = linspace(0, L, Nx);                  % discretize x axis domain
    [X,Y] = meshgrid(x,x) ;                 % create 2D grid of points
    Wmn = sin(m*pi*X/L) .* sin(n*pi*Y/L);    % mode shapes
    Wnm = sin(n*pi*X/L) .* sin(m*pi*Y/L);
end

% figure set-up
v     = [0,0] ;
figure1 = figure('Color',[0 0 0]);
colormap(bone); 

% plotting loop - animates the plate for varying theta
for th = 1:length(theta)
    % degenerate mode shape
    w = cos(theta(th))*Wmn + sin(theta(th))*Wnm;
    
    % plotting
    mesh(X,Y,w,'facealpha',0.3) ; hold on;
    set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',...
    [1 1 1]);
    contour(X,Y,w,v,'color','w','linewidth',2) ; hold off;
    zlim([-2 2])
    xlabel('x')
    ylabel('y')
    zlabel('w')
    drawnow;
end


