clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

ode_name = 'Duffing';

xrange = 0;
nse_id = -1; %% unnoised: -1

ode_name = 'Duffing';

eqn_name = [ode_name '_xrange' num2str(xrange)];

fam_name = 'Chebyshev1';

bor = 10;

eqn_name = [ode_name '_xrange' num2str(xrange)];

plt_spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

plt_spc.color = [0 0 0 0.2];
solspc_ref_plot = LD_plots.plot_solspc( Sref, ...
                            LD_plots('Sref',[5 5],[1 5],[1 1],1), ...
                            plt_spc);

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

Sref_svd_plot_array = Sref_svd_basic_array;

Sref_crv_svd_plot = LD_plots.plot_curve_svds(Sref_svd_plot_array, ...
                                LD_plots('svd',[9 9],[2 9],[9 1],1));
