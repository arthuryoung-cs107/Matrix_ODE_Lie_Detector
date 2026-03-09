clear
close all

plt0 = LD_plots('S', ...
                [3 4],...
                [1 2],[1 1], ...
                1);
plt0 = plt0.init_tiles_safe(1,2);
hold(plt0.axs, 'on');
box(plt0.axs,'on');
axs = plt0.axs;

%% simple Riccati
fdxu = @(x_,u_) u_./x_ - u_;
%% our Riccati
% fdxu = @(x_,u_) 2.0*u_./x_ - u_.*u_.*x_.*x_;

Fode = ode;
                        % casted as trivial vector field, e ~ x
Fode.ODEFcn = @(e_,s_) [ones(size(e_)) ; fdxu(s_(1,:), s_(2,:))];

x0 = 1e-1;
u0_vals = logspace(-1,1,10);
ef = 2.0; % epsilon varies from e0 = 0 to ef > 0
xf = x0 + ef;

% nevl = length(u0_vals);
nevl = 2*length(u0_vals);
epsevl = linspace(0.0,ef,nevl);
x_evl = epsevl+x0;
u_evl = nan(length(u0_vals),nevl );
for i = 1:length(u0_vals)
    Fode.InitialValue = [x0 ; u0_vals(i)];

    % phi_i = solve(Fode,0.0,ef);
    % phi_xu_i = phi_i.Solution;
    % phi_x_i = phi_xu_i(1,:);
    % phi_u_i = phi_xu_i(2,:);

    phi_i = solutionFcn(Fode,0.0,ef);
    phi_xu_i = phi_i(epsevl);
    phi_x_i = phi_xu_i(1,:);
    phi_u_i = phi_xu_i(2,:);

    u_evl(i,:) = phi_u_i;

    plot3(axs(1), phi_x_i, phi_u_i, fdxu(phi_x_i, phi_u_i), ...
                    'LineStyle', '-', ...
                    'LineWidth', 2, ...
                    'Marker', '.', ...
                    'MarkerSize', 10, ...
                    'Color', LD_plots.green4);
    plot(axs(2), phi_x_i, phi_u_i, ...
                    'LineStyle', '-', ...
                    'LineWidth', 1, ...
                    'Marker', '.', ...
                    'MarkerSize', 15, ...
                    'Color', LD_plots.green4);
end
set(axs, ...
    'YScale', 'linear',  ...
    'XScale', 'linear',  ...
    'TickLabelInterpreter','Latex', ...
    'FontSize',16 );
xlabel(axs, '$$ x $$', 'Interpreter','Latex','FontSize',16);
ylabel(axs, '$$ u $$', 'Interpreter','Latex','FontSize',16);
zlabel(axs(1), '$$ d_x u $$', 'Interpreter','Latex','FontSize',16);
view(axs(1), LD_plots.view_mat(6, :));

plt0.show_toolbar()

% pause

%% purely geometric surface representation
% [Xmesh,Umesh] = meshgrid( x_evl, u0_vals  );
[Xmesh,Umesh] = meshgrid( x_evl, linspace(min(u_evl(:)),max(u_evl(:)),nevl)  );
DXUmesh = reshape( fdxu( Xmesh(:),Umesh(:)) , size(Xmesh) );
surf(axs(1),Xmesh,Umesh,DXUmesh)
