clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

% Strue_obs = LD_observations_set(dir_name, [eqn_name '_true_obs'], dat_suff);

Strue_obs_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_obs'], dat_suff);
Strue_rec_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_rec'], dat_suff);
Strue_ext_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_ext'], dat_suff);

Strue_obs_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_obs'], dat_suff);
Strue_rec_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_rec'], dat_suff);
Strue_ext_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_ext'], dat_suff);

% Strue_ref = Strue_obs;
Strue_ref = Strue_obs_DoP853;

Sarray_full = [Strue_obs_DoPri5;Strue_ext_DoPri5;Strue_obs_DoP853;Strue_ext_DoP853];
Sarray = [Strue_obs_DoPri5;Strue_ext_DoPri5;Strue_ext_DoP853];


solspc_ref_plot = LD_plots.plot_n2q1_solspc(Strue_ref,LD_plots('Strue_ref',[4 4],[1 4],[1 1],1));
[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Strue_ref,LD_plots('rel_err',[4 4],[1 4],[2 1],1));
