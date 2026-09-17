%% Lorenz Attractor: Numerical Solution and Animation

clear;
close all;
% clc;

%% Lorenz parameters
sigma = 10;
rho   = 28;
beta  = 8/3;

%% Initial condition
x0 = [1; 1; 1];

%% Time interval
t0 = 0;
tf = 50;

% Output times.
% Specifying these explicitly gives us convenient, evenly spaced
% points for the animation.
dt = 0.01;
tspan = t0:dt:tf;

%% Lorenz equations
% y(1) = x
% y(2) = y
% y(3) = z

lorenz = @(t,y) [ ...
    sigma*(y(2) - y(1));
    y(1)*(rho - y(3)) - y(2);
    y(1)*y(2) - beta*y(3)
];

%% Solve the ODE
[t,Y] = ode45(lorenz, tspan, x0);

x = Y(:,1);
y = Y(:,2);
z = Y(:,3);

%% Set up figure
fig_trj = figure( ...
'Theme', 'dark', ...
'MenuBar', 'none' ...
);
% 'Color','w' ...

ax = axes;
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

xlabel(ax,'$ u_1 $', 'Interpreter','Latex','FontSize',20);
ylabel(ax,'$ u_2 $', 'Interpreter','Latex','FontSize',20);
zlabel(ax,'$ u_3 $', 'Interpreter','Latex','FontSize',20);

title(ax,'Lorenz Attractor');

view(ax,3);

% Set fixed axis limits so they do not change during animation
xlim(ax,[min(x)-2, max(x)+2]);
ylim(ax,[min(y)-2, max(y)+2]);
zlim(ax,[min(z)-2, max(z)+2]);

%% Animated trajectory
trajectory = animatedline( ...
    'Color', 0.5*[1 1 1], ...
    'LineWidth', 1.5 ...
    );

% Marker indicating the current position
particle = plot3( ...
    x(1), y(1), z(1), ...
    'o', ...
    'MarkerFaceColor',[1 1 1], ...
    'MarkerEdgeColor',[1 1 1], ...
    'MarkerSize',10);

%% Animation parameters

% Plot several numerical time steps per video frame.
% Increase this to make the animation faster.
% frameSkip = 3;
frameSkip = 4;
% frameSkip = 1;

title_i = title(ax, ...
sprintf('$ j = %d , x_j = %.2f $',1,0), ...
'Interpreter','Latex','FontSize',20 ...
);

% trjvid = VideoWriter([sim_name '.avi'], 'Motion JPEG AVI');
trjvid = VideoWriter([[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'] 'Lorenz' '.mp4'], 'MPEG-4');
open(trjvid)
j = 1;
for k = 1:frameSkip:length(t)

    % Add new points to trajectory
    i0 = max(1,k-frameSkip+1);

    addpoints(trajectory, ...
        x(i0:k), ...
        y(i0:k), ...
        z(i0:k));

    % Move particle
    set(particle, ...
        'XData',x(k), ...
        'YData',y(k), ...
        'ZData',z(k));

    % Display current simulation time
    % title(ax,sprintf('Lorenz Attractor   t = %.2f',t(k)));
    % title(ax,sprintf('t = %.2f',t(k)));
    title_i.String = sprintf('$ j = %d , x_j = %.2f $',j,t(k));
    j = j+1;

    writeVideo(trjvid,getframe(fig_trj));

    % pause(0.05)
    % drawnow
    % drawnow limitrate;

end
close(trjvid);
