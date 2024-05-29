clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

spc.color = [0 0 0 0.2];
solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

Sref_Lsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Lmat');
Sref_Rsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Rmat');
Sref_Osvd = Sref.read_mat_svd_package('Chebyshev1',10,'Omat');
Sref_Psvd = Sref.read_mat_svd_package('Chebyshev1',10,'Pmat');
Sref_Qsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Qmat');
Sref_Gsvd = Sref.read_mat_svd_package('Chebyshev1',10,'Gmat');
Sref_OGsvd = Sref.read_mat_svd_package('Chebyshev1',10,'OGmat');

Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Psvd, ...
                            Sref_Qsvd, ...
                            Sref_Gsvd, ...
                            Sref_OGsvd];

[Sref_RYLsvd,Sref_RYLsvd_rst] = Sref.read_mat_svd_package('Chebyshev1',10,'Rmat','YL','Lmat');
[Sref_QYLsvd,Sref_QYLsvd_rst] = Sref.read_mat_svd_package('Chebyshev1',10,'Qmat','YL','Lmat');

[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
[Sref_PYLsvd,Sref_PYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Psvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
[Sref_OGYLsvd,Sref_OGYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_OGsvd,Sref_Lsvd);

% Sref_svd_YL_array = Sref_RYLsvd;
% Sref_svd_YL_array = [Sref_RYLsvd, Sref_QYLsvd];
Sref_svd_YL_array = [       Sref_RYLsvd, ...
                            Sref_OYLsvd, ...
                            Sref_PYLsvd, ...
                            Sref_QYLsvd, ...
                            Sref_GYLsvd, ...
                            Sref_OGYLsvd];

Sref_svdarray = [Sref_svd_basic_array, Sref_svd_YL_array];

Sref_svd_plot = LD_plots.plot_S_svds(Sref_svdarray,LD_plots('svd',[7 7],[2 7],[5 1],1));

S1array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Rmat','DoP853','ext',dat_suff);
S1array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Qmat','DoP853','ext',dat_suff);

S2array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'RmatYL','DoP853','ext',dat_suff);
S2array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'QmatYL','DoP853','ext',dat_suff);

Sarray = [S1array S2array];

[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Sref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));
