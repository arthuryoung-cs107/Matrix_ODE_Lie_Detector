%{
    ----------------------------
    preliminaries
    ----------------------------
%}

clear;
close all;

sys_screens = apv_plots.get_sys_screens();
screen_specs = sys_screens(1,:);
[o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1)

% scrn_id = 1;
scrn_id = 3;

%{
    ----------------------------
    Define ODE system
    ----------------------------
%}

fode = struct( ...
'name', 'Riccati', ...
'eor', 1, ...
'ndep', 1, ...
'f', @(x_,u_) 2.0*( u_ ./ x_ ) - (x_.*x_).*(u_.*u_) , ...
'gradf', @(x_,u_) [ -2.0*(u_ ./ (x_.*x_)) - 2.0*(x_).*(u_.*u_) ; ...
2.0*(1.0./x_) - 2.0*(x_.*x_).*u_ ] ...
);

Fode_sys = ode;
Fode_sys.ODEFcn = @(e_,s_) [ ones(1,size(s_,2)) ; fode.f( s_(1,:),s_(2:end,:) ) ];

x0 = 1e-1;
u0_vals = logspace(-1,1,10);
ef = 2.0; % epsilon varies from e0 = 0 to ef > 0
xf = x0 + ef;

ncrv = length(u0_vals);
nevl = 33; % one more than the cubic Rmat curve matrix

epsevl = linspace(0.0,ef,nevl);
x_evl = epsevl+x0;
u_evl = nan(length(u0_vals),nevl );
Sobs = cell(ncrv,1);
for i = 1:length(u0_vals)
    Fode_sys.InitialValue = [x0 ; u0_vals(i)];

    phi_i = solutionFcn(Fode_sys,0.0,ef);
    phi_xu_i = phi_i(epsevl);
    u_evl(i,:) = phi_xu_i(2:end,:);

    Sobs{i} = [ phi_xu_i ; fode.f( phi_xu_i(1,:),phi_xu_i(2:end,:) ) ]; % sols are the graph of f
end

%{
    ----------------------------
    plot results of ode system
    ----------------------------
%}

plt0 = apv_plots('S', ...
                [3 4],...
                [1 2],[1 1], ...
                1);
plt0 = plt0.init_tiles_safe(1,2);
hold(plt0.axs, 'on');
box(plt0.axs,'on');
axs = plt0.axs;

plot3(axs(1), x_evl, u_evl, fode.f(x_evl, u_evl), ...
                'LineStyle', '-', ...
                'LineWidth', 2, ...
                'Marker', '.', ...
                'MarkerSize', 10, ...
                'Color', LD_plots.green4);
plot(axs(2), x_evl, u_evl, ...
                'LineStyle', '-', ...
                'LineWidth', 1, ...
                'Marker', '.', ...
                'MarkerSize', 15, ...
                'Color', LD_plots.green4);

set(axs, ...
    'YScale', 'linear',  ...
    'XScale', 'linear',  ...
    'TickLabelInterpreter','Latex', ...
    'FontSize',16 );
xlabel(axs, '$$ x $$', 'Interpreter','Latex','FontSize',16);
ylabel(axs, '$$ u $$', 'Interpreter','Latex','FontSize',16);
zlabel(axs(1), '$$ d_x u $$', 'Interpreter','Latex','FontSize',16);
view(axs(1), apv_plots.view_mat(6, :));




%{
    ----------------------------
    model observed ode system
    ----------------------------
%}

Fode_obs = fode; % debugging
% Fode_obs = struct( 'eor', 1, 'ndep', 1 );

% Fode_obs.get_s_icrv_jsol = @(S_,i_,j_) S_{i_}(:,j_);
Fode_obs.get_s_icrv_jsol = @(S_,i_,j_) S_{i_}(:,j_);
Fode_obs.unpack_Smat = @(s_) deal( s_(1,:),s_( 2:(Fode_obs.ndep+1),: ),s_( (Fode_obs.ndep+2):end,: ));
Fode_obs.unpack_Smat_xun = @(s_) deal( s_(1,:),s_( 2:(1+Fode_obs.ndep*(Fode_obs.eor+1)),: ) );
Fode_obs.unpack_s = @(s_) unpack_Smat( s_ );
Fode_obs.unpack_s_xun = @(s_) unpack_Smat_xun( s_ );

bor = 3;
% fspace0 = apv.make_polynomial_fspace(1+fode.ndep,bor)
fspace0 = apv(Fode_obs); % Fode_obs must have ndep, eor
fspace0 = fspace0.init_polynomial_fspace(bor) % initialize to mv polynomial

Smat_obs = jspc.Scell_2_Smat(Sobs,jspc.ndim(fspace0));

[Renc,specs] = apv.model_Fode_observations(fspace0,Fode_obs,Sobs)

fprintf('(end)\n');
