clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

Stru_8gen = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
Stru_5gen = LD_observations_set(dir_name,eqn_name,'true','DoPri5', dat_suff);
Sarray_gen = [Stru_8gen;Stru_5gen];
Strue_ref = Stru_8gen;
% Strue_ref = Stru_5gen;

solspc_ref_plot = LD_plots.plot_n2q1_solspc(Strue_ref,LD_plots('Strue_ref',[4 4],[1 4],[1 1],1));

SRarray(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Rmat','DoP853','ext',dat_suff);
SRarray(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',9,'Rmat','DoP853','ext',dat_suff);
SRarray(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',8,'Rmat','DoP853','ext',dat_suff);

SQarray(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Qmat','DoP853','ext',dat_suff);
SQarray(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',9,'Qmat','DoP853','ext',dat_suff);
SQarray(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',8,'Qmat','DoP853','ext',dat_suff);

Sarray = [SRarray SQarray];

S_svd_plot = LD_plots.plot_S_svds(Sarray,LD_plots('svd',[5 5],[1 5],[3 1],1));
[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Strue_ref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));
