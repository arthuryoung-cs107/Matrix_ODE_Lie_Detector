%{
    ----------------------------
    matlab preliminaries
    ----------------------------
%}
clear;
close all;

%{
    ----------------------------
    plotting contexts
    ----------------------------
%}
sys_screens = apv_plots.get_sys_screens();
screen_specs = sys_screens(1,:);
[o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1)
scrn_id = 1;

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}
% % generate a data set
tic0 = tic;
[Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Riccati_data_2(); dat_true.bor_max = 3; % N = 1, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Brusselator_data(); % N = 1, Q = 2
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Van_der_Pol_data(); % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_oscillator_polr_data(); % N = 2, Q = 1, easier than VanderPol
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_pendulum_polr_data(); % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_double_oscillator_data(); dat_true.bor_max = 3; % N = 2, Q = 2, linear homogenous
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_double_pendulum_data(); % N = 2, Q = 2, would be wild to learn anything
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Linden_bouyancy_data(); dat_true.bor_max = 6; % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Lorenz_data(); dat_true.bor_max = 3; % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Gilmore_data(); dat_true.bor_max = 3; % N = 2, Q = 1
toc1 = toc(tic0);
fprintf('generated jet space data in %.3f seconds \n', toc1);

ncrv = length(Sobs(:));
npts_per_crv = zeros(ncrv,1);
for i = 1:ncrv
    npts_per_crv(i) = size(Sobs{i},2);
end
icrv_check = 1;
isol_check = 10;
iisol_check = sum(npts_per_crv(1:(icrv_check-1))) + isol_check;

dat = struct('ndep',dat_true.ndep,'eor',dat_true.eor);
    %% choose an arbitrary origin for the generation of an intrinsic coordinate system
    dat.icrv_sO = 2;
        % dat.i_sO_0 = 1;
        dat.i_sO_0 = 2;
        % dat.i_sO_0 = ceil(npts_per_crv(dat.icrv_sO)/10);
        % dat.i_sO_0 = ceil(npts_per_crv(dat.icrv_sO)/5);
        % dat.i_sO_0 = ceil(npts_per_crv(dat.icrv_sO)/4);
        % dat.i_sO_0 = ceil(npts_per_crv(dat.icrv_sO)/3);
        % dat.i_sO_0 = ceil(npts_per_crv(dat.icrv_sO)/2);
        % dat.i_sO_0 = ceil(npts_per_crv(dat.icrv_sO)*(2/3));
    dat.i_sO_1 = dat.i_sO_0+1;

ndep = dat.ndep;
eor = dat.eor;
nvar = 1+ndep;
ndim = 1+ndep*(eor + 1);

%{
    ----------------------------
    plot observations
    ----------------------------
%}

plt0 = apv_plots('jetspace', ...
                [3 6],...
                [1 4],[2 2], ...
                scrn_id);
dat_plt0 = dat;
dat_plt0.LineStyle = '-';
dat_plt0.Color = apv_plots.green4;
dat_plt0.MarkerSize = 6;

plt0 = apv_plots.plot_Sobs(plt0,Sobs,dat_plt0);

dat_plt0.Color = apv_plots.blue1;
plt0 = apv_plots.plot_Sobs(plt0,Sobs{icrv_check},dat_plt0);
dat_plt0.Color = apv_plots.blue5;
plt0 = apv_plots.plot_Sobs(plt0,Sobs{icrv_check}(:,isol_check) ,dat_plt0);

dat_plt0.Color = apv_plots.green4;
plt0.show_toolbar

% keyboard
% return
%{
    ----------------------------
    specify function space
    ----------------------------
%}
Sall  = cell2mat(Sobs(:)');
lo = min(Sall(1:nvar,:),[],2);  hi = max(Sall(1:nvar,:),[],2);
fmap0 = struct('Omap_a', 2./(hi-lo), 'Omap_b', -(hi+lo)./(hi-lo));
% fmap0 = struct( ...
% 'Omap_a', ones(nvar,1), ...
% 'Omap_b', zeros(nvar,1) ...
% );
fmap = fmap0;

%{
    ----------------------------
    model solution space
    ----------------------------
%}
tic0 = tic;
mod = LDsol.model_solspace(Sobs,dat,fmap);
toc1 = toc(tic0);
fprintf('built jet space model in %.3f seconds \n', toc1);

mod
MarkerSize_0 = dat_plt0.MarkerSize;
dat_plt0.MarkerSize = ceil(0.8*MarkerSize_0);
% sO = mod.s_O
dat_plt0.Color = apv_plots.grey3;
plt0 = apv_plots.plot_Sobs(plt0, Sobs{mod.flow_pckg.icrv_nghbr},dat_plt0);
dat_plt0.Color = [1 1 1];
plt0 = apv_plots.plot_Sobs(plt0, Sobs{dat.icrv_sO},dat_plt0);
dat_plt0.MarkerSize = ceil(1.5*MarkerSize_0);

plt0 = apv_plots.plot_Sobs(plt0, mod.s_O ,dat_plt0);
dat_plt0.Color = apv_plots.grey2;
plt0 = apv_plots.plot_Sobs(plt0, mod.Smat(:,mod.flow_pckg.isol_nghbr) ,dat_plt0);
dat_plt0.Color = apv_plots.grey4;
plt0 = apv_plots.plot_Sobs(plt0, mod.Smat(:,mod.flow_pckg.isol_nghbr+1) ,dat_plt0);
dat_plt0.Color = apv_plots.orange3;
plt0 = apv_plots.plot_Sobs(plt0, mod.flow_pckg.s_E ,dat_plt0);
dat_plt0.Color = apv_plots.red5;
plt0 = apv_plots.plot_Sobs(plt0, mod.flow_pckg.s_D ,dat_plt0);

dat_plt0.Color = apv_plots.green4;
dat_plt0.MarkerSize = MarkerSize_0;

Jf_RN1_check = mod.J_tau_u_RN1(:,:,1,iisol_check)
Jdxf_RN1_check = mod.J_tau_u_RN1(:,:,2,iisol_check)
JF_tru_check = JF_obs{icrv_check}( :,:,isol_check)
JF_N1_check = mod.JF_N1( :,:, iisol_check )

dNp1xu_tru_check = dNp1xu_obs{icrv_check}( :,isol_check )
dNp1xu_N1mod_check = mod.tau_uN_RN1_net((end-ndep+1):end,2,iisol_check)

[JF_sO_tru,dNp1xu_sO_tru] = dat_true.JF_dxf(mod.s_O);

dNp1xu_sO_tru
dNp1xu_sO_check = mod.t_O((end-ndep+1):end)
err_dNp1xu_sO = abs((dNp1xu_sO_check-dNp1xu_sO_tru)./dNp1xu_sO_tru)

VN_sO_unit = mod.flow_pckg.VN_spc_unit_sO
VN1_sO_unit = mod.flow_pckg.VN1_spc_unit_sO
JF_Vspc_sO = mod.sol_O.JF * mod.flow_pckg.VN_spc_sO

norms_VN_sO = sqrt(sum(mod.flow_pckg.VN_spc_sO.^2,1))

UJFT_VN_sO = mod.flow_pckg.JF_sO_svd.U' * mod.flow_pckg.VN_spc_unit_sO
UJFT_VN1_sO = mod.flow_pckg.JF_N1_sO_svd.U' * mod.flow_pckg.VN1_spc_unit_sO
VNT_VN_sO = mod.flow_pckg.VN_spc_unit_sO' * mod.flow_pckg.VN_spc_unit_sO
VN1T_VN1_sO = mod.flow_pckg.VN1_spc_unit_sO' * mod.flow_pckg.VN1_spc_unit_sO

% mu_svd_sO = mod.flow_pckg.Gn_bse.mu_N1_sO_svd;
% check = mu_svd_sO.U * diag(mu_svd_sO.s) * mu_svd_sO.V';
% check(:,1:10)
Jxi_sO = mod.flow_pckg.Jxi_sO;
cond_Jxi_sO = cond(Jxi_sO)
rank_Jxi_sO = rank(Jxi_sO)
% V0T_Jxi_sO = mod.flow_pckg.VN_spc_unit_sO( 1:(1+ndep*eor),: )' * mod.flow_pckg.Jxi_sO
V0T_Jxi_sO = mod.flow_pckg.VN_spc_sO( 1:(1+ndep*eor),: )' * mod.flow_pckg.Jxi_sO

dat_plt1 = dat_plt0;
plt1 = apv_plots('model_summary', ...
                [1 1],...
                [1 1],[1 1], ...
                scrn_id);
% [1 1],...
% [1 1],[1 1], ...
% [plt1,plt2,plt3,plt4] = apv_plots.plot_LDsol_model_summary(plt1,mod,dat_plt1);
[plt1,plt2] = apv_plots.plot_LDsol_model_summary(plt1,mod,dat_plt1);
% plt1.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])
% plt2.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])
% plt3.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])
% plt4.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])
