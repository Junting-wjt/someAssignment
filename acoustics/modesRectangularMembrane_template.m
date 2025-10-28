%-----------------------------------------------------
%           Acoustics Tutorial 7 
%    Plotting the modes of a rectangular membrane (skeleton code)
%            Dr Thomas McKenzie
%           Reid School of Music 
%          University of Edinburgh   
%-----------------------------------------------------

clc
clear all
close all

% input params
Lx = 1;             % membrane length in x-direction
Ly = 1;             % membrane length in y-direction
m = 5; n = 4;       % modal index pair
bc = 0;             % choice of boundary condition -> 0: free; 1: fixed

% create grid
dx = 0.01*Lx;               % grid spacing in x
dy = 0.01*Ly;               % grid spacing in y
x = 0 : dx : Lx;            % x axis    
y = 0 : dy : Ly;            % y axis
[X,Y] = meshgrid(x,y) ;     % create 2D grid of points

% mode shapes...
Wmn = sin(m*pi*X/Lx) .* sin(n*pi*Y/Ly);
Wnm = sin(n*pi*X/Lx) .* sin(m*pi*Y/Ly);

% angle of node lines
theta = 2*pi*(0:99)/100;

% open figure
figure1 = figure('Color',[0 0 0]);


for th = 1:length(theta)
    
    % degenerate mode shape
    W = cos(theta(th))*Wmn + sin(theta(th))*Wnm;
    
    % plot
    v     = [0,0] ;
    colormap(bone); 
    mesh(X,Y,W,'facealpha',0.3) ; hold on;
    set(gca,'Color',[0 0 0],'XColor',[1 1 1],'YColor',[1 1 1],'ZColor',...
    [1 1 1]);
    contour(X,Y,W,v,'color','w','linewidth',2) ; hold off;
    zlim([-2 2])
    xlabel('x')
    ylabel('y')
    zlabel('w')
    drawnow
end

W = cos(theta).*Wmn + sin(theta).*Wnm;



