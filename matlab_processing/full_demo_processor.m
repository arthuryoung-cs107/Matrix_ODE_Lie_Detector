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

[~,s_Ofull,V_Ofull] = svd(Sref_Osvd.matT','econ','vector');

solspc_ref_plot = LD_plots.plot_n2q1_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1));
Sref_svd_plot = LD_plots.plot_S_svds(Sref_svdarray,LD_plots('svd',[5 5],[1 5],[3 1],1));

SRarray(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Rmat','DoP853','ext',dat_suff);
SRarray(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',9,'Rmat','DoP853','ext',dat_suff);
SRarray(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',8,'Rmat','DoP853','ext',dat_suff);

SQarray(1) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',10,'Qmat','DoP853','ext',dat_suff);
SQarray(2) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',9,'Qmat','DoP853','ext',dat_suff);
SQarray(3) = LD_observations_set(dir_name,eqn_name,'true','DoP853','Chebyshev1',8,'Qmat','DoP853','ext',dat_suff);

Sarray = [SRarray SQarray];

S_svd_plot = LD_plots.plot_S_svds(Sarray,LD_plots('svd',[5 5],[1 5],[4 1],1));
[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Sref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));