clear;
close all;

% dir_name = '../data_directory';
% dir_name = '../denoise_data_directory';
dir_name = '../denoise_data_directory/Gaussian_IC_perturbation';
dense_dir_name = '../dense_data_directory/Gaussian_IC_perturbation';
dat_suff = 'lddat';

xrange = 0;
% xrange = 1;

ode_name = 'Duffing';
% ode_name = 'VanDerPol';
% ode_name = 'Pendulum';
% ode_name = 'Bessel';
% ode_name = 'Riccati';
% ode_name = 'Brusselator';

eqn_name = [ode_name '_xrange' num2str(xrange)];
% eqn_name = [ode_name '_xrange' num2str(xrange) '_extrap'];

fam_name = 'Chebyshev1';
% fam_name = 'Chebyshev2';
% fam_name = 'Legendre';

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

% spc.color = [0 0 0 0.2];
% solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

bor = 10;
fspace0 = LD_orthopolynomial_space(bor,meta0);
fspace0 = fspace0.read_domain_configuration(Sref.make_fspace_config_name(fam_name,bor));

noise_level = 1;
Snse = LD_observations_set(dir_name,eqn_name,['noise' num2str(noise_level)],'DoP853', dat_suff);
% solspc_nse_plot = LD_plots.plot_solspc(Snse,LD_plots('Snse',[4 4],[1 4],[2 1],1),spc);

Sref.crvs = Sref.make_curve_array
Snse.crvs = Snse.make_curve_array

%% inspecting one curve
i_crv = 1;

tdim_S_mesh = 5;
tdim_S_hght = 3;
spc.color = [0 0 0 0.2];
nse_plot1 = LD_denoise_plots.plot_observed_trajectories(LD_plots('Snse', ...
                                                        [tdim_S_mesh tdim_S_mesh],[tdim_S_hght tdim_S_mesh],[tdim_S_hght 1],1), ...
                                                        spc, ...
                                                        Sref,Snse,i_crv);
% spc.color = [0 1 0 0.2];
% Sref2 = LD_observations_set('../data_directory',eqn_name,'true','DoP853', dat_suff);
% Snse2 = LD_observations_set('../data_directory',eqn_name,['noise' num2str(noise_level)],'DoP853', dat_suff);
% nse_plot2 = LD_denoise_plots.plot_observed_trajectories(nse_plot1, ...
%                                                         spc, ...
%                                                         Sref2,Snse2);

% rjet_ref = LD_observations_set.regularize_curve_jets(crvs_ref);


% A_name = 'L';
% A_name = 'G';
% A_name = 'Rn';
A_name = 'Rnp1';
% A_name = 'OG';


% L_ref = Sref.read_encoded_matrices([fam_name '.' num2str(bor) '.L']);
% G_ref = Sref.read_encoded_matrices([fam_name '.' num2str(bor) '.G']);
% Rn_ref = Sref.read_encoded_matrices([fam_name '.' num2str(bor) '.Rn']);
% Rnp1_ref = Sref.read_encoded_matrices([fam_name '.' num2str(bor) '.Rnp1']);
% OG_ref = Sref.read_encoded_matrices([fam_name '.' num2str(bor) '.OG']);

[A_ref,Asvd_ref] = Sref.read_encoded_matrices([fam_name '.' num2str(bor) '.' A_name]);
[A_nse,Asvd_nse] = Snse.read_encoded_matrices([fam_name '.' num2str(bor) '.' A_name]);;

% Sdns = LD_observations_set(dense_dir_name,eqn_name,'true','DoP853', dat_suff);
% crvs_dns = Sdns.make_curve_array
% crv_i_d = crvs_dns(i_crv)
% jet_i_d = crv_i_d.jets


%% plots

% plt_check = nse_plot1;
% axs_check = plt_check.axs;
% plot(axs_check(1),x_check,u_check,'- o','Color',[0 1 0],'MarkerSize',3,'LineWidth',2)

kor_u_hat = meta0.eor;

tdim_E_mesh = 7;
tdim_E_hght = 4;
spc.mspec = 'none';
spc.ms = 3;
spc.lspec = '-';
spc.lw = 1;
est_plot1 = LD_denoise_plots.plot_curve_estimates(nse_plot1, ...
    LD_plots('error_estimation', ...
    [tdim_E_mesh tdim_E_mesh],[tdim_E_hght tdim_E_mesh],[tdim_E_mesh 1],1), ...
    spc,Sref,Snse,i_crv);

% nse_plot1.show_menubar();
nse_plot1.show_toolbar();

% est_plot1.show_menubar();
est_plot1.show_toolbar();

% mat_plot1 = LD_denoise_plots.plot_Amat_data(LD_plots('Amat',[7 7],[2 7],[7 1],1), spc, ...
%     A_ref,Asvd_ref, ...
%     A_nse,Asvd_nse);

fprintf('      [END of denoise_demo_processor.m]\n')
