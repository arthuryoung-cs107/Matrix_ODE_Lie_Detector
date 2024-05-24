clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
Sref_Lsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Lmat');
Sref_Rsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Rmat');
Sref_Osvd = Sref.read_mat_svd_package('Chebyshev1',10,'Omat');
Sref_Psvd = Sref.read_mat_svd_package('Chebyshev1',10,'Pmat');
Sref_Qsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Qmat');
Sref_Gsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Gmat');

Sref_svdarray = [   Sref_Lsvd, ...
                    Sref_Rsvd, ...
                    Sref_Osvd, ...
                    Sref_Psvd, ...
                    Sref_Qsvd, ...
                    Sref_Gsvd];

solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1));
Sref_svd_plot = LD_plots.plot_S_svds(Sref_svdarray,LD_plots('svd',[7 7],[3 7],[5 1],1));

S1array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'Qmat','DoP853','ext',dat_suff);
% S1array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'Qmat','DoP853','ext',dat_suff);

S2array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Qmat','DoP853','ext',dat_suff);
S2array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev2',10,'Qmat','DoP853','ext',dat_suff);

Sarray = [S1array S2array];

[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Sref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));

solspc_div_plot = LD_plots.plot_S_icrv_divergence(Sref, Sref.read_mat_svd_package('Legendre',10,'Gmat'),2, ...
                                                        solspc_ref_plot,LD_plots('Sdiv',[5 5],[2 5],[5 1],1));
