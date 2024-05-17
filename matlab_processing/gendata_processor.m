clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

% Strue_obs = LD_observations_set(dir_name, [eqn_name '_true_obs'], dat_suff);
% Strue_rec = LD_observations_set(dir_name, [eqn_name '_true_rec'], dat_suff);
% Strue_ext = LD_observations_set(dir_name, [eqn_name '_true_ext'], dat_suff);

Strue_obs_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_obs'], dat_suff);
Strue_rec_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_rec'], dat_suff);
Strue_ext_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_ext'], dat_suff);

Strue_obs_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_obs'], dat_suff);
Strue_rec_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_rec'], dat_suff);
Strue_ext_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_ext'], dat_suff);

% solspc_obs_plot = LD_plots.plot_n2q1_solspc(Strue_obs);
% solspc_rec_plot = LD_plots.plot_n1q1_solspc(Strue_rec);
% solspc_ext_plot = LD_plots.plot_n2q1_solspc(Strue_ext);

solspc_obs_DoPri5_plot = LD_plots.plot_n2q1_solspc(Strue_obs_DoPri5);
solspc_rec_DoPri5_plot = LD_plots.plot_n1q1_solspc(Strue_rec_DoPri5);
solspc_ext_DoPri5_plot = LD_plots.plot_n2q1_solspc(Strue_ext_DoPri5);

solspc_obs_DoP853_plot = LD_plots.plot_n2q1_solspc(Strue_obs_DoP853);
solspc_rec_DoP853_plot = LD_plots.plot_n1q1_solspc(Strue_rec_DoP853);
solspc_ext_DoP853_plot = LD_plots.plot_n2q1_solspc(Strue_ext_DoP853);
