clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

Sref_Lsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Lmat');
Sref_Rsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Rmat');
Sref_Osvd = Sref.read_mat_svd_package('Chebyshev1',10,'Omat');
Sref_Psvd = Sref.read_mat_svd_package('Chebyshev1',10,'Pmat');
Sref_Qsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Qmat');
Sref_Gsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Gmat');

Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Psvd, ...
                            Sref_Qsvd, ...
                            Sref_Gsvd];

Sref_RYLsvd = Sref.read_mat_svd_package('Chebyshev1',10,'RmatYL');
Sref_QYLsvd = Sref.read_mat_svd_package('Chebyshev1',10,'QmatYL');

Sref_svd_YL_array = [Sref_RYLsvd, Sref_QYLsvd];

Sref_svdarray = [Sref_svd_basic_array, Sref_svd_YL_array];

Sref_svd_plot = LD_plots.plot_S_svds(Sref_svdarray,LD_plots('svd',[7 7],[2 7],[5 1],1));

S1array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Rmat','DoP853','ext',dat_suff);
S1array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'RmatYL','DoP853','ext',dat_suff);
S1array(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Qmat','DoP853','ext',dat_suff);
S1array(4) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'QmatYL','DoP853','ext',dat_suff);

S2array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'Rmat','DoP853','ext',dat_suff);
S2array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'RmatYL','DoP853','ext',dat_suff);
S2array(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'Qmat','DoP853','ext',dat_suff);
S2array(4) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'QmatYL','DoP853','ext',dat_suff);

Sarray = [S1array S2array];

[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Sref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));

return

solspc_div_plot = LD_plots.plot_S_icrv_divergence(Sref, Sref.read_mat_svd_package('Chebyshev1',10,'Rmat'),1, ...
                                                        solspc_ref_plot,LD_plots('Sdiv',[5 5],[2 5],[5 1],1));
