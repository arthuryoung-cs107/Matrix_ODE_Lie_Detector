clear;
close all;

% eqn_name = 'Duffing';
eqn_name = 'Duffing_DoPri5';
dir_name = '../data_directory';
dat_suff = 'lddat';

Strue = LD_observations_set(dir_name, [eqn_name '_true_obs'], dat_suff);
Strue_rec = LD_observations_set(dir_name, [eqn_name '_true_rec'], dat_suff);
Strue_ext = LD_observations_set(dir_name, [eqn_name '_true_ext'], dat_suff);

solspc_plot = LD_plots.plot_n2q1_solspc(Strue);
solspc_rec_plot = LD_plots.plot_n1q1_solspc(Strue_rec);
solspc_ext_plot = LD_plots.plot_n2q1_solspc(Strue_ext);
