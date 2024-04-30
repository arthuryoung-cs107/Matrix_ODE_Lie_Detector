clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

Strue = LD_observations_set(dir_name, 'Duffing_true_obs', dat_suff);

solspc_plot = LD_plots.plot_n2q1_solspc(Strue);
