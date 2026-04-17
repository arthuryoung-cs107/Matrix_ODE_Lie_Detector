%{
    ----------------------------
    preliminaries
    ----------------------------
%}

clear;
close all;

sys_screens = apv_plots.get_sys_screens();
% sys_screens = get(groot,'MonitorPositions');

screen_specs = sys_screens(1,:);
[o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1)

scrn_id = 1;
% scrn_id = 3;

%{
    ----------------------------
    Define ODE system
    ----------------------------
%}
eqn_name = 'Riccati'
f_eqn = @(x_,u_) 2.0*( u_ ./ x_ ) - (x_.*x_).*(u_.*u_) ;
gradf_eqn = @(x_,u_) [ ...
    -2.0*(u_ ./ (x_.*x_)) - 2.0*(x_).*(u_.*u_) ; ...
    2.0*(1.0./x_) - 2.0*(x_.*x_).*u_ ...
];
dxf_eqn = @(x_,u_) sum( ...
    gradf_eqn( x_, u_(1,:) )' ...
    .* [ ones(size(x_))' , u_(2,:)' ] ...
    , 2  )';

fode = struct( ...
'name', eqn_name, ...
'eor', 1, ...
'ndep', 1, ...
'f', @(x_,u_) f_eqn(x_,u_), ...
'gradf', @(x_,u_) gradf_eqn(x_,u_), ...
'dxf', @(x_,u_) dxf_eqn(x_,u_) ...
);

xu_check = [ 1e-1 ; 1e1 ]
xu_ad = adobj(xu_check,eye(length(xu_check)))
% x_ad = xu_ad(1) % too glitchy, not worth it
% u_ad = xu_ad(2) % too glitchy, not worth it
x_ad = xu_ad.qdim(1)
u_ad = xu_ad.qdim(2)

% xu_ad_check = xu_ad(:)
% xu_ad_check = xu_ad(1:2)

f_check = fode.f(xu_check(1),xu_check(2))
g_check = fode.gradf(xu_check(1),xu_check(2))

f_ad = 2.0.*( u_ad ./ x_ad ) - (x_ad.*x_ad).*(u_ad.*u_ad)
fmap_ad = @(x_,u_) 2.0*( u_ ./ x_ ) - (x_.*x_).*(u_.*u_);

% fm_ad = fmap_ad( xu_ad(1),xu_ad(2) ) % too glitchy, not worth it
fm_ad = fmap_ad( xu_ad.qdim(1),xu_ad.qdim(2) )

if ( fm_ad.Jac(:) ~= g_check(:) )
    fprintf('autodiff is broken\n');
    pause
end

% Fode_sys = ode;
% Fode_sys.ODEFcn = @(e_,s_) [ ones(1,size(s_,2)) ; fode.f( s_(1,:),s_(2:end,:) ) ];
Fode_sys_evl = @(e_,s_) [ ones(1,size(s_,2)) ; f_eqn( s_(1,:),s_(2:end,:) ) ];

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
    % Fode_sys = ode;
    % Fode_sys.ODEFcn =

    Fode_sys = ode45(Fode_sys_evl,[0.0,ef],[x0 ; u0_vals(i)]);

    % phi_i = solutionFcn(Fode_sys,0.0,ef);
    % phi_xu_i = phi_i(epsevl);
    phi_xu_i = deval(Fode_sys,epsevl);

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
                [1 3],[1 1], ...
                1);
plt0 = plt0.init_tiles_safe(1,3);
hold(plt0.axs, 'on');
box(plt0.axs,'on');
axs = plt0.axs;

plot(axs(1), x_evl, u_evl, ...
    'LineStyle', '-', ...
    'LineWidth', 0.5, ...
    'Color', (apv_plots.green4), ...
    'Marker', 'o', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', apv_plots.green4, ...
    'MarkerEdgeColor', [0 0 0] ...
    );
plot(axs(2), x_evl, fode.f(x_evl, u_evl), ...
    'LineStyle', '-', ...
    'LineWidth', 0.5, ...
    'Color', (apv_plots.green4), ...
    'Marker', 'o', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', apv_plots.green4, ...
    'MarkerEdgeColor', [0 0 0] ...
    );
plot3(axs(3), x_evl, u_evl, fode.f(x_evl, u_evl), ...
    'LineStyle', '-', ...
    'LineWidth', 0.5, ...
    'Color', (apv_plots.green4), ...
    'Marker', 'o', ...
    'MarkerSize', 6, ...
    'MarkerFaceColor', apv_plots.green4, ...
    'MarkerEdgeColor', [0 0 0] ...
    );

set(axs, ...
    'YScale', 'linear',  ...
    'XScale', 'linear',  ...
    'TickLabelInterpreter','Latex', ...
    'FontSize',16 );
xlabel(axs, '$$ x $$', 'Interpreter','Latex','FontSize',16);
ylabel(axs, '$$ u $$', 'Interpreter','Latex','FontSize',16);
zlabel(axs(end), '$$ d_x u $$', 'Interpreter','Latex','FontSize',16);
view(axs(end), apv_plots.view_mat(6, :));




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

N1cubic = ldaux.first_order_ld_ad_cubic_model(Sobs,Fode_obs)

pause

% [tvf,Renc,specs] = apv.model_trivial_vfield( fspace0,Fode_obs,Sobs )
[mod,tvf,Renc,specs] = apv.model_Fode_observations( fspace0,Fode_obs,Sobs );

% mod
% tvf
% Renc
% specs

obj_ = fspace0;

fprintf('\nmod.Pbor_cell{1} * mod.G1Pbor_svd(1).V(:,end)')
fspc.print_vshort_polynomial_theta( mod.Pbor_cell{1} * mod.G1Pbor_svd(1).V(:,end) , obj_.P_mat);
fprintf('mod.Pbor_cell{1} * mod.G1Pbor_svd(1).W(:,end)')
fspc.print_vshort_polynomial_theta( mod.Pbor_cell{1} * mod.G1Pbor_svd(1).W(:,end) , obj_.P_mat);

fprintf('\nmod.Pbor_cell{2} * mod.G1Pbor_svd(2).V(:,end)')
fspc.print_vshort_polynomial_theta( mod.Pbor_cell{2} * mod.G1Pbor_svd(2).V(:,end) , obj_.P_mat);
fprintf('mod.Pbor_cell{2} * mod.G1Pbor_svd(2).W(:,end)')
fspc.print_vshort_polynomial_theta( mod.Pbor_cell{2} * mod.G1Pbor_svd(2).W(:,end) , obj_.P_mat);

% pre = '\nmod.W(:,end)';
% fspc.print_vshort_polynomial_theta( mod.W(:,end) , obj_.P_mat);
% fprintf('mod.W(:,end-1)')
% fspc.print_vshort_polynomial_theta( mod.W(:,end-1) , obj_.P_mat);
% fprintf('mod.W(:,end-2)')
% fspc.print_vshort_polynomial_theta( mod.W(:,end-2) , obj_.P_mat);
fprintf('\nmod.V(:,end)')
fspc.print_vshort_polynomial_theta( mod.V(:,end) , obj_.P_mat);
fprintf('mod.V(:,end-1)')
fspc.print_vshort_polynomial_theta( mod.V(:,end-1) , obj_.P_mat);
fprintf('mod.V(:,end-2)')
fspc.print_vshort_polynomial_theta( mod.V(:,end-2) , obj_.P_mat);

fprintf('\nmod.KGnet(:,1)')
fspc.print_vshort_polynomial_theta( mod.KGnet(:,1) , obj_.P_mat);
fprintf('mod.KGnet(:,2)')
fspc.print_vshort_polynomial_theta( mod.KGnet(:,2) , obj_.P_mat);

return

fprintf('\nmod.R1G1_svd.W(:,end)')
fspc.print_vshort_polynomial_theta(mod.R1G1_svd.W(:,end),obj_.P_mat);
fprintf('mod.R1G1_svd.W(:,end-1)')
fspc.print_vshort_polynomial_theta(mod.R1G1_svd.W(:,end-1),obj_.P_mat);

fprintf('\nmod.KGmKR_svd.U(:,1)')
fspc.print_vshort_polynomial_theta(mod.KGmKR_svd.U(:,1),obj_.P_mat);
fprintf('mod.WGmWR_svd.U(:,1)')
fspc.print_vshort_polynomial_theta(mod.WGmWR_svd.U(:,1),obj_.P_mat);

return

obj_ = fspace0;
            fspc.print_vshort_polynomial_theta(R1_svd.W(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(R1_svd.W(:,end-1),obj_.P_mat);


            % tvf.W_R1
            [R1_svd.Wx,R1_svd.Wu] = obj_.split_Vxu_mat(R1_svd.W);
            WR1x_svd = fspc.compute_svd_package(R1_svd.Wx)

            %{
                Scaled singular vectors of the Kx ~ Wx submatrix.
                Contains Kx subspace generating vartheta image
            %}
            fspc.print_vshort_polynomial_theta_z(WR1x_svd.Y(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta_z(WR1x_svd.Y(:,2),obj_.P_mat);

            %{
                Scaled singular vectors of the Vx submatrix.
            %}
            fspc.print_vshort_polynomial_theta(vartheta_svd.Y(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(vartheta_svd.Y(:,2),obj_.P_mat);

            rstats = @(r_) deal( sqrt(sum(r_)), max(r_(:)), median(r_(:)) );

            [tx_S,tu_S] = deal(txu_S(1,:), txu_S(2:end,:));
            res_tu_S = sum((tu_S-d1xu).^2,1);
            [err_tu_S_net,maxr_tu_S,medr_tu_S] = rstats(res_tu_S)

            K_R1 = vartheta_svd.V(:,1:vartheta_svd.r); % an orthonormal basis for ker(R1)
            [ txu_K_S , vartheta_K_S ] = obj_.txu(obj_,xu,K_R1); % pass back to txu

            [tx_K_S,tu_K_S] = deal(txu_K_S(1,:), txu_K_S(2:end,:));
            res_tu_K_S = sum((tu_K_S-d1xu).^2,1);
            [err_tu_K_S_net,maxr_tu_K_S,medr_tu_K_S] = rstats(res_tu_K_S)

            tu_mat_S = reshape(tu_S,ndep,nobs);
            vartheta_tns_S = reshape(vartheta_K_S,Plen,nvar,nobs);

            tdxu_mat_S = nan(ndep,nobs); % estimates of second derivative of u via tvf's first jetspace coord
            T0_full_mat = nan(ntheta,nobs);
            [T1_full_mat,H1_full_mat] = deal(nan(ntheta,nobs));
            G1_full_ttns = zeros(ntheta,ndep,nobs);
            G1_true_ttns = zeros(ntheta,ndep,nobs);
            for i = 1:nobs
                vtx_i = vartheta_tns_S(:,1,i);
                vtu_i = vartheta_tns_S(:,2:end,i);

                l_i = l_S(i,:)';
                g_i = gxu_l_S(:,:,i)';
                Lam1_x_i = reshape(Lam_x_dxu(:,:,1,i),ndep,Plen)';
                Lam1_u_i = Lam_u_dxu(:,1,i);

                Lam_dxu_T_i = G1_full_ttns(:,:,i);
                Lam0_xu_T_i = zeros(ntheta,nvar);

                % tdxu_mat_S = ;
                % T1_full_mat(1:Plen,:,i) = l_i;

                T0_full_mat(1:Plen,i) = l_i; % tx*vx = 1*vx
                Lam0_xu_T_i(1:Plen,1) = l_i;
                idel_P = Plen;
                for idep = 1:ndep
                    Lam_dxu_T_i(1:Plen,idep) = Lam1_x_i(:,idep);

                    inds_ii = (1+idel_P):(Plen+idel_P);
                    Lam_dxu_T_i(inds_ii,idep) = Lam1_u_i;
                    Lam0_xu_T_i(inds_ii,idep+1) = l_i;
                    T0_full_mat(inds_ii,i) = d1xu(:,idep)*l_i; % tuq*vuq = dxuq*li
                    idel_P = idel_P + Plen;
                end

                G1_full_ttns(:,:,i) = ...
                ((vtu_i' - (tu_mat_S(:,i).*vtx_i'))*g_i*Lam0_xu_T_i' - Lam_dxu_T_i')';

                % oracle verification
                G1_true_ttns(:,:,i) = ...
                (obj_.obs.gradf(Smat(1,i),Smat(2:nvar,i))'*Lam0_xu_T_i' - Lam_dxu_T_i')';
            end

            T0_svd = fspc.compute_svd_package(T0_full_mat)

            fspc.print_vshort_polynomial_theta(T0_svd.V(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(T0_svd.V(:,end-1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(T0_svd.V(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(T0_svd.V(:,2),obj_.P_mat);

            G1_svd = fspc.compute_svd_package(G1_full_ttns)

            fspc.print_vshort_polynomial_theta(G1_svd.W(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(G1_svd.W(:,end-1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(G1_svd.W(:,end-2),obj_.P_mat);

            G1_true_svd = fspc.compute_svd_package(G1_true_ttns)

            check = sqrt(sum( (G1_true_svd.mat*G1_svd.V(:,(end-3):end)).^2 ,1))

            fspc.print_vshort_polynomial_theta(G1_true_svd.W(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(G1_true_svd.W(:,end-1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(G1_true_svd.W(:,end-2),obj_.P_mat);

            %{
                The G matrix SVD reveals WG ~ KG, the kernel of G.
                The kernel of G forms a basis for all Lambda vector fields
                on the solution space, S

                It can immediately be used to refine the tvf model via projection
                onto vector field coordinates.
            %}

            WR_T_WG = ( R1_svd.W )' * (G1_svd.W);
            WR_T_WG_svd = fspc.compute_svd_package(WR_T_WG)

            KR_T_KG = ( R1_svd.PK(R1_svd) )' * (G1_svd.PK(G1_svd));
            KR_T_KG_svd = fspc.compute_svd_package(KR_T_KG)

            D_KR_T_KG = KR_T_KG_svd.D(KR_T_KG_svd);
            F_RG = G1_svd.PK(G1_svd)*D_KR_T_KG;

            fspc.print_vshort_polynomial_theta(F_RG(:,1),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(F_RG(:,2),obj_.P_mat);

            YRYG_svd = fspc.compute_svd_package( [ (G1_svd.Y)' ; (R1_svd.Y)' ] )
            fspc.print_vshort_polynomial_theta(YRYG_svd.V(:,end),obj_.P_mat);
            fspc.print_vshort_polynomial_theta(YRYG_svd.V(:,end-1),obj_.P_mat);



            % size(s_T1_full)
            % YT1_full = V_T1_full.*(s_T1_full');
            % YG1_full = V_G1_full.*(s_G1_full');
            % size(YT1_full)
            % size(YG1_full)
            % [W_T1G1,rank_T1G1,s_T1G1,V_T1G1,T1G1_full_mat] = ...
            %     fspc.safely_process_net_svd( [ YT1_full, YG1_full ] );
            % size(T1G1_full_mat)
            % rank_T1G1
            % s_T1G1_row = s_T1G1'
            % fspc.print_vshort_polynomial_theta(W_T1G1(:,end),obj_.P_mat);
            % fspc.print_vshort_polynomial_theta(W_T1G1(:,end-1),obj_.P_mat);

fprintf('(end)\n');
