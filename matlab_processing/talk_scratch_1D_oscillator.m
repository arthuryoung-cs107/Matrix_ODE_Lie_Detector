%% 1D Harmonic Oscillator: Numerical Solution and Animation

clear;
close all;

%% Oscillator parameter
omega = 2;

%% Initial condition
% y(1) = position
% y(2) = velocity
y0 = [1; 0];

%% Time interval
t0 = 0;
tf = 20;

dt = 0.001;
tspan = t0:dt:tf;

%% Harmonic oscillator equations
oscillator = @(t,y) [ ...
    y(2);
    -omega^2*y(1)
];

%% Solve the ODE
[t,Y] = ode45(oscillator, tspan, y0);

x = Y(:,1);
v = Y(:,2);

%% Set up figure
fig_trj = figure( ...
'Theme', 'dark', ...
'MenuBar', 'none' ...
);


ax = axes;
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

xlabel(ax, '$ u $', 'Interpreter','Latex','FontSize',16);
title(ax,'1D Harmonic Oscillator');

% Give the oscillator some room
xmax = 1.2*max(abs(x));

xlim(ax,[-xmax xmax]);
ylim(ax,[-1 1]);

% Draw equilibrium position
plot(ax,[0 0],[-0.2 0.2],'k--');

% Draw horizontal line
plot(ax,[-xmax xmax],[0 0],'k-');

%% Oscillating particle
particle = plot( ...
    ax,x(1),0, ...
    'o', ...
    'MarkerFaceColor',[1 1 1], ...
    'MarkerEdgeColor',[1 1 1], ...
    'MarkerSize',12);

%% Animation
frameSkip = 1;
% frameSkip = 2;

% trjvid = VideoWriter([sim_name '.avi'], 'Motion JPEG AVI');
% trjvid = VideoWriter([[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'] 'Lorenz' '.mp4'], 'MPEG-4');
% open(trjvid)
for k = 1:frameSkip:length(t)

    % Move particle
    set(particle, ...
        'XData',x(k), ...
        'YData',0);

    % Display current simulation time
    title(ax,sprintf('t = %.2f',t(k)));

    pause(0.001)

    drawnow limitrate;
end
% close(trjvid);
