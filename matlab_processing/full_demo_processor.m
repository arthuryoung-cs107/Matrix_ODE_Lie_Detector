clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

% xrange = 0;
xrange = 1;
ode_name = 'Duffing';

eqn_name = [ode_name '_xrange' num2str(xrange)];
% eqn_name = [ode_name '_xrange' num2str(xrange) '_extrap'];

% fam_name = 'Chebyshev1';
fam_name = 'Legendre';
bor = 10;

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

spc.color = [0 0 0 0.2];
solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

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


[Sref_RYLsvd_in,Sref_RYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Rmat','YL','Lmat');
% [Sref_QYLsvd_in,Sref_QYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Qmat','YL','Lmat');
% Sref_svd_YL_in_array = [    Sref_RYLsvd_in,Sref_QYLsvd_in];
Sref_svd_YL_in_array = [    Sref_RYLsvd_in];

[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
[Sref_PYLsvd,Sref_PYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Psvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
[Sref_OGYLsvd,Sref_OGYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_OGsvd,Sref_Lsvd);
Sref_svd_YL_array = [       Sref_OYLsvd, Sref_PYLsvd, Sref_GYLsvd,Sref_OGYLsvd];

Sref_OregRsvd = Sref.read_mat_svd_package(fam_name,bor,'OregRmat');
Sref_svd_Oreg_array = [ Sref_OregRsvd   ];

Sref_svd_array = [Sref_svd_basic_array, Sref_svd_YL_in_array, Sref_svd_YL_array, Sref_svd_Oreg_array];
Sref_crv_svd_plot = LD_plots.plot_curve_svds(Sref_svd_array,LD_plots('svd',[7 7],[2 7],[6 1],1));

S1array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Rmat','DoP853','ext',dat_suff);
S1array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'Rmat','DoP853','ext',dat_suff);
S1array(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',9,'Rmat','DoP853','ext',dat_suff);
S1array(4) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',9,'Rmat','DoP853','ext',dat_suff);

S2array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'RmatYL','DoP853','ext',dat_suff);
S2array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'RmatYL','DoP853','ext',dat_suff);
S2array(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',9,'RmatYL','DoP853','ext',dat_suff);
S2array(4) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',9,'RmatYL','DoP853','ext',dat_suff);

S3array(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'OregRmat','DoP853','ext',dat_suff);
% S3array(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Legendre',10,'OregRmatYL','DoP853','ext',dat_suff);

% Sarray = [S1array S2array];
Sarray = [S1array S2array S3array];

[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Sref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));
% err_tns_check = reshape(sum(S_relerr_results.abs_rel_diff,1),[],Sref.ncrv,size(S_relerr_results.abs_rel_diff,3));
err_tns_check = reshape(sum(S_relerr_results.abs_rel_diff(2:(meta0.ndep+1),:,:),1),[],Sref.ncrv,size(S_relerr_results.abs_rel_diff,3));
crverr_check = reshape(sum(err_tns_check,1),Sref.ncrv,size(S_relerr_results.abs_rel_diff,3));
[~,isort_check] = sort(sum(crverr_check,2));

% icrv_check = 1
% icrv_check = isort_check(1)
icrv_check = isort_check(end)

[~,iS_check] = max(crverr_check(icrv_check,:))
% [~,iS_check] = min(crverr_check(icrv_check,:));

spc.color = [0 0 0];
spc.lw = 2.0;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_check);

spc.color = [1 0 0];
spc.lw = 2.0;
solspc_ref_plot = LD_plots.plot_solspc(Sarray(iS_check),solspc_ref_plot,spc,icrv_check);
