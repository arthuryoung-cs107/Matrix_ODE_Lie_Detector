clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

nse_name = 'true';
gen_name = 'DoP853';
fam_name = 'Chebyshev1';
bor = 10;
mat_name = 'Rmat';
rec_name = 'DoP853';

asmbl_gen_name = @(nse_,gen_) [eqn_name '_' nse_ '_' gen_ 'gen'];
asmbl_syn_name = @(nse_,gen_,fam_,bor_,mat_,rec_) [asmbl_gen_name(nse_,gen_) '_' fam_ '.' num2str(bor_) '.' mat_ '.' rec_];

tru_8gen_name = asmbl_gen_name(nse_name,gen_name);
tru_8gen_Cheb1R10_8rec_name = asmbl_syn_name(nse_name,gen_name,fam_name,bor,mat_name,rec_name)

% Strue_obs = LD_observations_set(dir_name, [eqn_name '_true_obs'], dat_suff);

Strue_obs_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_obs'], dat_suff);
Strue_rec_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_rec'], dat_suff);
Strue_ext_DoPri5 = LD_observations_set(dir_name, [eqn_name '_DoPri5_true_ext'], dat_suff);

Strue_obs_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_obs'], dat_suff);
Strue_rec_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_rec'], dat_suff);
Strue_ext_DoP853 = LD_observations_set(dir_name, [eqn_name '_DoP853_true_ext'], dat_suff);

Stru_8gen = LD_observations_set(dir_name, tru_8gen_name, dat_suff);
Stru_8gen_Cheb1R10_8rec = LD_observations_set(dir_name, [tru_8gen_Cheb1R10_8rec_name 'rec'], dat_suff);
Stru_8gen_Cheb1R10_8ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R10_8rec_name 'ext'], dat_suff);

% Strue_ref = Strue_obs;
Strue_ref = Strue_obs_DoP853;

Sarray_full = [Strue_obs_DoPri5;Strue_ext_DoPri5;Strue_obs_DoP853;Strue_ext_DoP853;Stru_8gen;Stru_8gen_Cheb1R10_8ext];
Sarray = [Strue_obs_DoPri5;Strue_ext_DoPri5;Strue_ext_DoP853;Stru_8gen;Stru_8gen_Cheb1R10_8ext];


solspc_ref_plot = LD_plots.plot_n2q1_solspc(Strue_ref,LD_plots('Strue_ref',[4 4],[1 4],[1 1],1));
solspc_check_plot = LD_plots.plot_n2q1_solspc(Stru_8gen,LD_plots('S_check',[4 4],[1 4],[3 1],1));
[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Strue_ref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));
