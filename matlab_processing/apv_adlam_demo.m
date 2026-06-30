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
% sys_screens = get(groot,'MonitorPositions');

screen_specs = sys_screens(1,:);
[o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1)

scrn_id = 1;
% scrn_id = 2;
% scrn_id = 3;

%{
    ----------------------------
    load/generate observations
    ----------------------------
%}
% load a data set


% generate a data set
tic0 = tic;
% [Sobs,dat_true] = ldaux.generate_Riccati_data(); % N = 1, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Riccati_data_2(); % N = 1, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Brusselator_data(); % N = 1, Q = 2
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_Van_der_Pol_data(); % N = 2, Q = 1
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_oscillator_polr_data(); % N = 2, Q = 1, easier than VanderPol
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_oscillator_cart_data(); % N = 2, Q = 2 via map
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_pendulum_polr_data(); % N = 2, Q = 1
% [Sobs,dat_true] = ldaux.generate_pendulum_cart_data(); % N = 2, Q = 2
[Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_double_oscillator_data(); % N = 2, Q = 2, linear homogenous
% [Sobs,dat_true,JF_obs,dNp1xu_obs] = ldaux.generate_double_pendulum_data(); % N = 2, Q = 2, would be wild to learn anything
toc1 = toc(tic0);
fprintf('generated jet space data in %.3f seconds \n', toc1);
s0 = Sobs{1}(:,1);

% JF_obs
% dNp1xu_obs

icrv_check = 1;
isol_check = 10;

s_check = Sobs{icrv_check}(:,isol_check);
sNm1 = s_check(1:(end-dat_true.ndep));
JF_tru_check = JF_obs{icrv_check}(:,:,isol_check)

% f_sNm1 = dat_true.f(sNm1);
% Jf_sNm1 = dat_true.Jf(sNm1);
% s_check
% sNm1
% f_sNm1
% Jf_sNm1
% sNm1_ad = adfcn.seed_sol(s_check(1:(end-dat_true.ndep)));
sNm1_ad = adobj.seed_sol(sNm1);
% sNm1_ad
% f_ad_sNm1 = adobj.f_ad( sNm1 , dat_true.f_ad );
% f_ad_sNm1
% f_ad_sNm1.Jac

% return
% f_ad_sNm1 = adfcn.f_ad( s_check( 1:(end-dat_true.ndep+1) ) , dat_true.f_ad );

% dat = dat_true; % oracle information, debugging purposes
dat = struct('ndep',dat_true.ndep,'eor',dat_true.eor);

ndep = dat.ndep;
eor = dat.eor;
nvar = 1+ndep;
ndim = 1+ndep*(eor + 1)
%{
    ----------------------------
    plot observations
    ----------------------------
%}
dat_plt0 = dat;
dat_plt0.LineStyle = '-';
dat_plt0.Color = apv_plots.green4;
plt0 = apv_plots('jetspace', ...
                [3 6],...
                [1 4],[2 2], ...
                scrn_id);
% [3 4],...
% [1 3],[1 1], ...
% [2 1],...
% [1 1],[1 1], ...
% plt0 = apv_plots.plot_Sobs(plt0,Sobs,dat_plt0);
plt0 = apv_plots.plot_Sobs(plt0,Sobs{icrv_check},dat_plt0);
dat_plt0.Color = apv_plots.blue1;
plt0 = apv_plots.plot_Sobs(plt0,Sobs{icrv_check}(:,isol_check) ,dat_plt0);
dat_plt0.Color = apv_plots.green4;

%{
    ----------------------------
    build jet space model
    ----------------------------
%}
fspace0 = struct( ...
'bor', 3, ...
'Omap_a', ones(1+dat.ndep,1), ...
'Omap_b', zeros(1+dat.ndep,1) ...
);
fspace = fspace0;
% fspace.Omap_b = -s0(1:(dat.ndep+1)); % set origin to first observed solution, w.l.o.g.
% fspace.fam = 'Legendre';
% fspace.fam = 'Hermite1';
% fspace.fam = 'Hermite2';
% fspace.fam = 'Chebyshev1';
% fspace.fam = 'Chebyshev2';
% fspace = mvp_jspc_model.Omap_solutions(fspace,Sobs,[-1,1]);
% fspace = mvp_jspc_model.Omap_solutions(fspace,Sobs,[0,1]);
% fspace = struct( ...
% 'bor', 3, ...
% 'Omap_a', [ 0.1 ; [ 1.1 ; 1.2 ].*ones(dat.ndep,1) ], ...
% 'Omap_b', zeros(1+dat.ndep,1) + [1.0 ; 2.0; 3.0] ...
% ); % debugging
tic0 = tic;
mod = mvp_jspc_model(Sobs,dat,fspace)
toc1 = toc(tic0);
fprintf('built jet space model in %.3f seconds \n', toc1);

%{
    ----------------------------
    review jet space model
    ----------------------------
%}
mod = mvp_jspc_model.verify(mod);

iisol_check = sum(mod.npts_per_crv(1:(icrv_check-1))) + isol_check

Jf_N1mod_check = mod.jspc_N1mod.J_tau_u_glb(:,:,1,iisol_check)
Jdxf_N1mod_check = mod.jspc_N1mod.J_tau_u_glb(:,:,2,iisol_check)
JF_tru_check = JF_obs{icrv_check}( :,:,isol_check)
JF_N1mod_check = mod.jspc_N1mod.JF_glb( :,:, iisol_check )

dNp1xu_tru_check = dNp1xu_obs{icrv_check}( :,isol_check )
dNp1xu_N1mod_check = mod.jspc_N1mod.tau_uN_net((end-ndep+1):end,2,iisol_check)

%{
    ----------------------------
    plot jet space model data
    ----------------------------
%}
% dat_plt1.LineStyle = '-';
% dat_plt1.Color = apv_plots.green4;
plt1 = apv_plots('model_summary', ...
                [1 1],...
                [1 1],[1 1], ...
                scrn_id);
plt1 = apv_plots.plot_model_summary(plt1,mod,dat_plt0);

return
