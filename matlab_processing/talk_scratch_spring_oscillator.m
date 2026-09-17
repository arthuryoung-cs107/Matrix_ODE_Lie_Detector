%% 1D Mass-Spring Oscillator Animation

clear;
close all;

%% Oscillator parameters
omega = 2;
c_damp = 2e-1;

%% Initial conditions
% y(1) = displacement
% y(2) = velocity
y0 = [0.5; 0];

%% Time interval
t0 = 0;
tf = 20;

% dt = 0.001;
dt = 0.01;
tspan = t0:dt:tf;

%% Harmonic oscillator
oscillator = @(t,y) [ ...
    y(2);
    -omega^2*y(1)-c_damp*y(2)
];

%% Solve
[t,Y] = ode45(oscillator,tspan,y0);

x = Y(:,1);
v = Y(:,2);

%% Geometry of animation

% Equilibrium vertical position of mass
y_eq = -1.5;

% Convert oscillator displacement into vertical position
y_mass = y_eq + x;

% Mass dimensions
massWidth  = 0.35;
massHeight = 0.25;

% Suspension point
y_support = 0;

%% Set up figure
fig_trj = figure( ...
'Theme', 'dark', ...
'MenuBar', 'none' ...
);

ax = axes;
hold(ax,'on');
grid(ax,'on');
box(ax,'on');

xlim(ax,[-1 1]);
ylim(ax,[-2.5 0.5]);

% xlabel(ax,'$ $');
ylabel(ax,'$ u $', 'Interpreter','Latex','FontSize',20);

title(ax,'1D Mass-Spring Oscillator');

ax.XTick = [];

%% Draw fixed support

supportWidth  = 0.6;
supportHeight = 0.12;

rectangle( ...
    'Position', ...
    [-supportWidth/2, ...
      y_support, ...
      supportWidth, ...
      supportHeight], ...
    'FaceColor',[0.7 0.7 0.7], ...
    'EdgeColor','k', ...
    'LineWidth',1.5);

%% Draw dashed suspension line

lineHandle = plot( ...
    [0 0], ...
    [y_support y_mass(1)+massHeight/2], ...
    'Color', [1 1 1], ...
    'LineStyle', '-.', ...
    'LineWidth',2);

%% Draw mass

mass = rectangle( ...
    'Position',[ ...
        -massWidth/2, ...
        y_mass(1)-massHeight/2, ...
        massWidth, ...
        massHeight], ...
    'FaceColor',[0.1 0.4 0.9], ...
    'EdgeColor','k', ...
    'LineWidth',1.5);

%% Animation

frameSkip = 4;

title_i = title(ax, ...
sprintf('$ j = %d , x_j = %.2f $',1,0), ...
'Interpreter','Latex','FontSize',20 ...
);

trjvid = VideoWriter([[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'] 'spring_mass' '.mp4'], 'MPEG-4');
open(trjvid)
j = 1;
for k = 1:frameSkip:length(t)

    % Current mass center
    yc = y_mass(k);

    % Update mass
    mass.Position = [ ...
        -massWidth/2, ...
        yc-massHeight/2, ...
        massWidth, ...
        massHeight];

    % Update suspension line
    lineHandle.YData = [ ...
        y_support, ...
        yc+massHeight/2];

    % Update title
    % title(ax,sprintf( ...
    %     ' t = %.2f', ...
    %     t(k)));
    title_i.String = sprintf('$ j = %d , x_j = %.2f $',j,t(k));
    j = j+1;

    % drawnow limitrate;

    writeVideo(trjvid,getframe(fig_trj));
end
close(trjvid);
