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

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}
% % generate a data set
tic0 = tic;
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Riccati_data_2(); % N = 1, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Brusselator_data(); % N = 1, Q = 2
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Van_der_Pol_data(); % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_oscillator_polr_data(); % N = 2, Q = 1, easier than VanderPol
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_pendulum_polr_data(); % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_double_oscillator_data(); % N = 2, Q = 2, linear homogenous
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_double_pendulum_data(); % N = 2, Q = 2, would be wild to learn anything
[Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Linden_bouyancy_data(); % N = 2, Q = 2, linear homogenous
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
plt0 = apv_plots.plot_Sobs(plt0,Sobs,dat_plt0);
dat_plt0.Color = apv_plots.blue1;
plt0 = apv_plots.plot_Sobs(plt0,Sobs{icrv_check},dat_plt0);
dat_plt0.Color = apv_plots.blue5;
plt0 = apv_plots.plot_Sobs(plt0,Sobs{icrv_check}(:,isol_check) ,dat_plt0);
dat_plt0.Color = apv_plots.green4;
plt0.show_toolbar

% return
%{
    ----------------------------
    specify function space
    ----------------------------
%}
fmap0 = struct( ...
'Omap_a', ones(nvar,1), ...
'Omap_b', zeros(nvar,1) ...
);
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
sO = mod.s_O
dat_plt0.Color = apv_plots.red1;
plt0 = apv_plots.plot_Sobs(plt0, Sobs{mod.isrtmags_dXi_S_sO(2)},dat_plt0);
dat_plt0.Color = [1 1 1];
plt0 = apv_plots.plot_Sobs(plt0, sO ,dat_plt0);
dat_plt0.Color = apv_plots.green4;

Jf_RN1_check = mod.J_tau_u_RN1(:,:,1,iisol_check)
Jdxf_RN1_check = mod.J_tau_u_RN1(:,:,2,iisol_check)
JF_tru_check = JF_obs{icrv_check}( :,:,isol_check)
JF_N1_check = mod.JF_N1( :,:, iisol_check )

dNp1xu_tru_check = dNp1xu_obs{icrv_check}( :,isol_check )
dNp1xu_N1mod_check = mod.tau_uN_RN1_net((end-ndep+1):end,2,iisol_check)

dat_plt1 = dat_plt0;
plt1 = apv_plots('model_summary', ...
                [1 1],...
                [1 1],[1 1], ...
                scrn_id);
% [1 1],...
% [1 1],[1 1], ...
plt1 = apv_plots.plot_LDsol_model_summary(plt1,mod,dat_plt1);
