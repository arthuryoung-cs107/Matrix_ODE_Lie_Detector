clear;
close all;

eqn_name = 'Duffing';
dir_name = '../data_directory';
dat_suff = 'lddat';

asmbl_gen_name = @(nse_,gen_) [eqn_name '_' nse_ '_' gen_ 'gen'];
asmbl_syn_name = @(nse_,gen_,fam_,bor_,mat_,rec_) [asmbl_gen_name(nse_,gen_) '_' fam_ '.' num2str(bor_) '.' mat_ '.' rec_];

tru_8gen_name = asmbl_gen_name('true','DoP853')
tru_8gen_Cheb1R10_8rec_name = asmbl_syn_name('true','DoP853','Chebyshev1',10,'Rmat','DoP853')
tru_8gen_Cheb1R10_5rec_name = asmbl_syn_name('true','DoP853','Chebyshev1',10,'Rmat','DoPri5')
tru_8gen_Cheb1R8_8rec_name = asmbl_syn_name('true','DoP853','Chebyshev1',8,'Rmat','DoP853')
tru_8gen_Cheb1R8_5rec_name = asmbl_syn_name('true','DoP853','Chebyshev1',8,'Rmat','DoPri5')

tru_5gen_name = asmbl_gen_name('true','DoPri5')
tru_5gen_Cheb1R10_8rec_name = asmbl_syn_name('true','DoPri5','Chebyshev1',10,'Rmat','DoP853')
tru_5gen_Cheb1R10_5rec_name = asmbl_syn_name('true','DoPri5','Chebyshev1',10,'Rmat','DoPri5')
tru_5gen_Cheb1R8_8rec_name = asmbl_syn_name('true','DoPri5','Chebyshev1',8,'Rmat','DoP853')
tru_5gen_Cheb1R8_5rec_name = asmbl_syn_name('true','DoPri5','Chebyshev1',8,'Rmat','DoPri5')
tru_5gen_Cheb1R10_8rec_name = asmbl_syn_name('true','DoPri5','Legendre',10,'Rmat','DoP853')
tru_5gen_Cheb1R10_5rec_name = asmbl_syn_name('true','DoPri5','Legendre',10,'Rmat','DoPri5')
tru_5gen_Cheb1R8_8rec_name = asmbl_syn_name('true','DoPri5','Legendre',8,'Rmat','DoP853')
tru_5gen_Cheb1R8_5rec_name = asmbl_syn_name('true','DoPri5','Legendre',8,'Rmat','DoPri5')

Stru_8gen = LD_observations_set(dir_name, tru_8gen_name, dat_suff);
Stru_8gen_Cheb1R10_8ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R10_8rec_name 'ext'], dat_suff);
Stru_8gen_Cheb1R10_5ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R10_5rec_name 'ext'], dat_suff);
Stru_8gen_Cheb1R8_8ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R8_8rec_name 'ext'], dat_suff);
Stru_8gen_Cheb1R8_5ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R8_5rec_name 'ext'], dat_suff);
Stru_8gen_LegenR10_8ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R10_8rec_name 'ext'], dat_suff);
Stru_8gen_LegenR10_5ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R10_5rec_name 'ext'], dat_suff);
Stru_8gen_LegenR8_8ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R8_8rec_name 'ext'], dat_suff);
Stru_8gen_LegenR8_5ext = LD_observations_set(dir_name, [tru_8gen_Cheb1R8_5rec_name 'ext'], dat_suff);

Stru_5gen = LD_observations_set(dir_name, tru_5gen_name, dat_suff);
Stru_5gen_Cheb1R10_8ext = LD_observations_set(dir_name, [tru_5gen_Cheb1R10_8rec_name 'ext'], dat_suff);
Stru_5gen_Cheb1R10_5ext = LD_observations_set(dir_name, [tru_5gen_Cheb1R10_5rec_name 'ext'], dat_suff);
Stru_5gen_Cheb1R8_8ext = LD_observations_set(dir_name, [tru_5gen_Cheb1R8_8rec_name 'ext'], dat_suff);
Stru_5gen_Cheb1R8_5ext = LD_observations_set(dir_name, [tru_5gen_Cheb1R8_5rec_name 'ext'], dat_suff);
Stru_5gen_LegenR10_8ext = LD_observations_set(dir_name, [tru_8gen_LegenR10_8rec_name 'ext'], dat_suff);
Stru_5gen_LegenR10_5ext = LD_observations_set(dir_name, [tru_8gen_LegenR10_5rec_name 'ext'], dat_suff);
Stru_5gen_LegenR8_8ext = LD_observations_set(dir_name, [tru_8gen_LegenR8_8rec_name 'ext'], dat_suff);
Stru_5gen_LegenR8_5ext = LD_observations_set(dir_name, [tru_8gen_LegenR8_5rec_name 'ext'], dat_suff);

Strue_ref = Stru_8gen;
Sarray = [  Stru_8gen; ...
            Stru_8gen_Cheb1R10_8ext; ...
            Stru_8gen_Cheb1R10_5ext; ...
            Stru_8gen_Cheb1R8_8ext; ...
            Stru_8gen_Cheb1R8_5ext; ...
            Stru_8gen_LegenR10_8ext; ...
            Stru_8gen_LegenR10_5ext; ...
            Stru_8gen_LegenR8_8ext; ...
            Stru_8gen_LegenR8_5ext; ...
            Stru_5gen; ...
            Stru_5gen_Cheb1R10_8ext; ...
            Stru_5gen_Cheb1R10_5ext; ...
            Stru_5gen_Cheb1R8_8ext; ...
            Stru_5gen_Cheb1R8_5ext; ...
            Stru_5gen_LegenR10_8ext; ...
            Stru_5gen_LegenR10_5ext; ...
            Stru_5gen_LegenR8_8ext; ...
            Stru_5gen_LegenR8_5ext];

solspc_ref_plot = LD_plots.plot_n2q1_solspc(Strue_ref,LD_plots('Strue_ref',[4 4],[1 4],[1 1],1));
[S_relerr_plot,S_relerr_results] = LD_plots.plot_S_relative_errors(Sarray,Strue_ref,LD_plots('rel_err',[5 5],[1 5],[5 1],1));
