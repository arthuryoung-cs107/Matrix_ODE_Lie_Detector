clear;
close all;

% dir_name = '../data_directory';
dir_name = '../denoise_data_directory';
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

bor = 10;

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

spc.color = [0 0 0 0.2];
solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);
fspace0 = LD_orthopolynomial_space(bor,meta0);
fspace0 = fspace0.read_domain_configuration(Sref.make_fspace_config_name(fam_name,bor));

noise_level = 2;
Snse = LD_observations_set(dir_name,eqn_name,['noise' num2str(noise_level)],'DoP853', dat_suff);

solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);
