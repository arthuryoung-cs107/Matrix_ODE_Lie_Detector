clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

fam_name = 'Legendre';
% fam_name = 'Chebyshev1';
% fam_name = 'Chebyshev2';

bor = 10;

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);
Sref_Lsvd = Sref.read_mat_svd_package(fam_name,bor,'Lmat');
Sref_Rsvd = Sref.read_mat_svd_package(fam_name,bor,'Rmat');
Sref_Osvd = Sref.read_mat_svd_package(fam_name,bor,'Omat');
Sref_Psvd = Sref.read_mat_svd_package(fam_name,bor,'Pmat');
Sref_Qsvd = Sref.read_mat_svd_package(fam_name,bor,'Qmat');
Sref_Gsvd = Sref.read_mat_svd_package(fam_name,bor,'Gmat');

Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Psvd, ...
                            Sref_Qsvd, ...
                            Sref_Gsvd];

[Sref_RYLsvd_in,Sref_RYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Rmat','YL','Lmat');
[Sref_QYLsvd_in,Sref_QYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Qmat','YL','Lmat');
Sref_svd_YL_in_array = [Sref_RYLsvd_in, Sref_QYLsvd_in];

% [Sref_RYLsvd,Sref_RYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Rsvd,Sref_Lsvd);
[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
% [Sref_PYLsvd,Sref_PYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Psvd,Sref_Lsvd);
% [Sref_QYLsvd,Sref_QYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Qsvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
% Sref_svd_YL_array = [       Sref_RYLsvd, ...
%                             Sref_OYLsvd, ...
%                             Sref_PYLsvd, ...
%                             Sref_QYLsvd, ...
%                             Sref_GYLsvd];
Sref_svd_YL_array = [       Sref_OYLsvd, ...
                            Sref_GYLsvd];

% Sref_svdarray = [Sref_svd_basic_array, Sref_svd_YL_in_array, Sref_svd_YL_array];
Sref_svdarray = [Sref_svd_basic_array, Sref_svd_YL_in_array, Sref_svd_YL_array];

Sref_svd_plot = LD_plots.plot_S_svds_w_global(Sref_svdarray,LD_plots('svd',[7 7],[3 7],[7 1],1));



% solspc_div_plot = LD_plots.plot_S_icrv_divergence(Sref, Sref.read_mat_svd_package('Chebyshev1',10,'Qmat'),1, ...
%                                                         solspc_ref_plot,LD_plots('Sdiv',[5 5],[2 5],[5 1],1));
% null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref.read_mat_svd_package('Chebyshev1',10,'Omat'),1, ...
%                                                         solspc_ref_plot,LD_plots('Sdiv',[5 5],[2 5],[5 1],1));
