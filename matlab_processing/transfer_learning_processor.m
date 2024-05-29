clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

% fam_name = 'Legendre';
fam_name = 'Chebyshev1';
% fam_name = 'Chebyshev2';

bor = 10;

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

spc.color = [0 0 0 0.2];
solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

fspace0 = LD_orthopolynomial_space(bor,meta0);
fspace0 = fspace0.read_domain_configuration(Sref.make_fspace_config_name(fam_name,bor));

Sref_Lsvd = Sref.read_mat_svd_package(fam_name,bor,'Lmat');
Sref_Rsvd = Sref.read_mat_svd_package(fam_name,bor,'Rmat');
Sref_Osvd = Sref.read_mat_svd_package(fam_name,bor,'Omat');
% Sref_Psvd = Sref.read_mat_svd_package(fam_name,bor,'Pmat');
% Sref_Qsvd = Sref.read_mat_svd_package(fam_name,bor,'Qmat');
Sref_Gsvd = Sref.read_mat_svd_package(fam_name,bor,'Gmat');
Sref_OGsvd = Sref.read_mat_svd_package(fam_name,bor,'OGmat');

Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Gsvd, ...
                            Sref_OGsvd];

[Sref_RYLsvd_in,Sref_RYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Rmat','YL','Lmat');
Sref_svd_YL_in_array = Sref_RYLsvd_in;

[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
[Sref_OGYLsvd,Sref_OGYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_OGsvd,Sref_Lsvd);
Sref_svd_YL_array = [       Sref_OYLsvd, ...
                            Sref_GYLsvd, ...
                            Sref_OGYLsvd];

Sref_svdarray = [Sref_svd_basic_array, Sref_svd_YL_in_array, Sref_svd_YL_array];

Sref_svd_plot = LD_plots.plot_S_svds_w_global(Sref_svdarray,LD_plots('svd',[7 7],[3 7],[7 1],1));

% solspc_div_plot = LD_plots.plot_S_icrv_divergence(Sref, Sref_Gsvd,2, ...
%                                                         solspc_ref_plot,LD_plots('Sdiv',[5 5],[2 5],[5 1],1));
% return
icrv_check = 4;
null_crv_plot = LD_plots('Snull',[5 5],[2 5],[5 1],1);
null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref_Gsvd,icrv_check, ...
                                                        solspc_ref_plot,null_crv_plot, LD_plots.get_green_mat);
null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref_Osvd,icrv_check, ...
                                                        solspc_ref_plot,null_crv_plot, LD_plots.get_blue_mat);
% null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref_OGsvd,icrv_check, ...
%                                                         solspc_ref_plot,null_crv_plot, LD_plots.get_purple_mat);
