clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

% xrange = 0;
xrange = 1;
ode_name = 'Duffing';

eqn_name = [ode_name '_xrange' num2str(xrange)];
% eqn_name = [ode_name '_xrange' num2str(xrange) '_extrap'];

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
Sref_Psvd = Sref.read_mat_svd_package(fam_name,bor,'Pmat');
Sref_Qsvd = Sref.read_mat_svd_package(fam_name,bor,'Qmat');
Sref_Gsvd = Sref.read_mat_svd_package(fam_name,bor,'Gmat');
Sref_OPsvd = Sref.read_mat_svd_package(fam_name,bor,'OPmat');
Sref_OGsvd = Sref.read_mat_svd_package(fam_name,bor,'OGmat');
Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Psvd, ...
                            Sref_Qsvd, ...
                            Sref_Gsvd, ...
                            Sref_OPsvd, ...
                            Sref_OGsvd];

Sref_svd_glb_basic_array = LD_aux.make_global_svd_package(Sref_svd_basic_array);
Sref_Lsvd_glb = Sref_svd_glb_basic_array(1);
Sref_Rsvd_glb = Sref_svd_glb_basic_array(2);
Sref_Osvd_glb = Sref_svd_glb_basic_array(3);
Sref_Psvd_glb = Sref_svd_glb_basic_array(4);
Sref_Qsvd_glb = Sref_svd_glb_basic_array(5);
Sref_Gsvd_glb = Sref_svd_glb_basic_array(6);
Sref_OPsvd_glb = Sref_svd_glb_basic_array(7);
Sref_OGsvd_glb = Sref_svd_glb_basic_array(8);


Sref_glb_svd_plot = LD_plots.plot_global_svds(Sref_svd_glb_basic_array,LD_plots('svd',[7 7],[2 7],[4 1],1));
Sref_glb_svd_plot.show_toolbar();


[Sref_RYLsvd_in,Sref_RYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Rmat','YL','Lmat');
% [Sref_QYLsvd_in,Sref_QYLsvd_rst_in] = Sref.read_mat_svd_package('Chebyshev1',10,'Qmat','YL','Lmat');
% Sref_svd_YL_in_array = [    Sref_RYLsvd_in,Sref_QYLsvd_in];
Sref_svd_YL_in_array = [    Sref_RYLsvd_in];

[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
[Sref_PYLsvd,Sref_PYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Psvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
[Sref_OGYLsvd,Sref_OGYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_OGsvd,Sref_Lsvd);
Sref_svd_YL_array = [       Sref_OYLsvd, Sref_PYLsvd, Sref_GYLsvd,Sref_OGYLsvd];

Sref_OregRsvd = Sref.read_mat_svd_package('Chebyshev1',10,'OregRmat');
Sref_svd_Oreg_array = [ Sref_OregRsvd   ];

Sref_svd_array = [Sref_svd_basic_array, Sref_svd_YL_in_array, Sref_svd_YL_array, Sref_svd_Oreg_array];
Sref_crv_svd_plot = LD_plots.plot_curve_svds(Sref_svdarray,LD_plots('svd',[7 7],[2 7],[6 1],1));

[~,isort_rvec] = sort(Sref_Gsvd.rvec);
% icrv_check = 2;
% icrv_check = isort_rvec(1);
icrv_check = isort_rvec(25);
% icrv_check = isort_rvec(end-22);
% solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_Osvd_glb,Sref_Lsvd_glb,icrv_check,fspace0, ...
% solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_OYLsvd_rst,Sref_Lsvd_glb,icrv_check,fspace0, ...
% solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_Rsvd_glb,Sref_Lsvd_glb,icrv_check,fspace0, ...
% solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_QYLsvd_rst_in,Sref_Lsvd_glb,icrv_check,fspace0, ...
% solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_Qsvd_glb,Sref_Lsvd_glb,icrv_check,fspace0, ...
% solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_Rsvd_glb,Sref_Lsvd_glb,icrv_check,fspace0, ...
solspc_glb_plot = LD_plots.plot_global_component(Sref,  Sref_Osvd_glb,Sref_RYLsvd_rst_in,Sref_Lsvd_glb,icrv_check,fspace0, ...
                                                        solspc_ref_plot,LD_plots('Sglb',[7 7],[1 7],[7 1],1));
solspc_glb_plot.show_toolbar();

% solspc_div_plot = LD_plots.plot_S_icrv_divergence(Sref, Sref_Gsvd,2, ...
%                                                         solspc_ref_plot,LD_plots('Sdiv',[5 5],[2 5],[5 1],1));

% null_crv_plot = LD_plots('Snull',[5 5],[2 5],[5 1],1);
% null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref_Gsvd,icrv_check, ...
%                                                         solspc_ref_plot,null_crv_plot, LD_plots.get_green_mat);
% null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref_Osvd,icrv_check, ...
%                                                         solspc_ref_plot,null_crv_plot, LD_plots.get_blue_mat);
% null_crv_plot = LD_plots.plot_nullity_icrv_comp(Sref,   Sref_OGsvd,icrv_check, ...
%                                                         solspc_ref_plot,null_crv_plot, LD_plots.get_purple_mat);
