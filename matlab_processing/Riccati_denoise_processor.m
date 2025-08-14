clear;
close all;

scrn_id = 1;

% dir_name = '../data_directory';
% dir_name = '../denoise_data_directory';
% dir_name = '../denoise_data_directory/Gaussian_IC_perturbation';
dir_name = '../denoise_data_directory/Gaussian_IC_perturbation/rendering_data';
% dense_dir_name = '../dense_data_directory/Gaussian_IC_perturbation';
dat_suff = 'lddat';

xrange = 0;
% xrange = 1;

ode_name = 'Riccati';

eqn_name = [ode_name '_xrange' num2str(xrange)];
% eqn_name = [ode_name '_xrange' num2str(xrange) '_extrap'];

% fam_name = 'Chebyshev1';
% fam_name = 'Chebyshev2';
fam_name = 'Legendre';

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

% spc.color = [0 0 0 0.2];
% solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

bor = 3;
fspace0 = LD_orthopolynomial_space(bor,meta0);
fspace0 = fspace0.read_domain_configuration(Sref.make_fspace_config_name(fam_name,bor));

% noise_level = 0;
noise_level = 1;
% noise_level = 2;
Snse = LD_observations_set(dir_name,eqn_name,['noise' num2str(noise_level)],'DoP853', dat_suff);
% solspc_nse_plot = LD_plots.plot_solspc(Snse,LD_plots('Snse',[4 4],[1 4],[2 1],1),spc);

Sref.crvs = Sref.make_curve_array;
Snse.crvs = Snse.make_curve_array;

Snse1 = LD_observations_set(dir_name,eqn_name,['noise' num2str(1)],'DoP853', dat_suff);
% Snse2 = LD_observations_set(dir_name,eqn_name,['noise' num2str(2)],'DoP853', dat_suff);

Lsvd_gf1jet = Sref.read_LD_svd('Lsvd_global_f1jet');
R1svd_gf1jet = Sref.read_LD_svd('R1svd_global_f1jet');
    Lsvd_gf1jet_h = Sref.read_LD_svd('Lsvd_global_f1jet_h');
    R1svd_gf1jet_h = Sref.read_LD_svd('R1svd_global_f1jet_h');
Rnsvd_g = Sref.read_LD_svd('Rnsvd_global',fam_name,bor);
rwimg_Rn_g = Sref.read_rowspace_image('Rnsvd_global',fam_name,bor);

Lsvd_gf1jet_n1 = Snse1.read_LD_svd('Lsvd_global_f1jet');
R1svd_gf1jet_n1 = Snse1.read_LD_svd('R1svd_global_f1jet');
    Lsvd_gf1jet_n1h = Snse1.read_LD_svd('Lsvd_global_f1jet_h');
    R1svd_gf1jet_n1h = Snse1.read_LD_svd('R1svd_global_f1jet_h');
Rnsvd_g_n1 = Snse1.read_LD_svd('Rnsvd_global',fam_name,bor);
rwimg_Rn_g_n1 = Snse1.read_rowspace_image('Rnsvd_global',fam_name,bor);
rwimg_Rn_st_g_n1 = Snse1.read_rowspace_image('Rnsvd_strue_global',fam_name,bor);
%
% Lsvd_gf1jet_n2 = Snse2.read_LD_svd('Lsvd_global_f1jet');
% R1svd_gf1jet_n2 = Snse2.read_LD_svd('R1svd_global_f1jet');
%     Lsvd_gf1jet_n2h = Snse2.read_LD_svd('Lsvd_global_f1jet_h');
%     R1svd_gf1jet_n2h = Snse2.read_LD_svd('R1svd_global_f1jet_h');
% Rnsvd_g_n2 = Snse2.read_LD_svd('Rnsvd_global',fam_name,bor);
% rwimg_Rn_g_n2 = Snse2.read_rowspace_image('Rnsvd_global',fam_name,bor);
% rwimg_Rn_st_g_n2 = Snse2.read_rowspace_image('Rnsvd_strue_global',fam_name,bor);

%% inspecting one curve
i_crv = 1;
% i_crv = 2;
% i_crv = 3;

% Snse_dns1 = LD_observations_set(dir_name,eqn_name,['noise' num2str(noise_level)],'DoP853','.jsol_R1_1', dat_suff);

tdim_S_mesh = 5;
tdim_S_wdth = 5;
tdim_S_hght = 2;
spc.color = [0 0 0 0.2];
% nse_plot1 = LD_denoise_plots.plot_observed_trajectories(LD_plots('Snse', ...
nse_plot1 = LD_denoise_plots.plot_denoised_trajectories(LD_plots('Sdns', ...
                                                        [tdim_S_mesh tdim_S_mesh],[tdim_S_hght tdim_S_wdth],[tdim_S_hght 1],scrn_id), ...
                                                        spc, ...
                                                        Sref,Snse,i_crv);
return
jet_sol_names = { '.jsol_h'; ...
'.jsol_h_R1'; ...
'.jsol_0_R1'; ...
'.jsol_1_R1'; ...
};
% '.jsol_R1' ; ...
[theta_jsh,pSjh,pSjhR1,pSj0R1,pSj1R1] = Snse.read_jet_sol_h_data('.theta_mat', jet_sol_names);
% [theta_jsh,pSjh,pSjhR1,pSj0R1,pSj1R1,pSjR1] = Snse.read_jet_sol_h_data('.theta_mat', jet_sol_names);
% [theta_jsh_1,pSjh_1,pSjhR1_1,pSj0R1_1,pSj1R1_1,pSjR1_1] = Snse_dns1.read_jet_sol_h_data('.theta_mat', jet_sol_names);

get_pts_mat = @(pS_) reshape(pS_.pts_in,meta0.ndim,[]);
pts_mat_crvi = @(pS_,i_) reshape(pS_.pts_in( (pS_.pts_crv_inds(1,i_)):(pS_.pts_crv_inds(2,i_)) ) , meta0.ndim,[]);

[spc.mspec,spc.ms] = deal('.',5);
[spc.lspec,spc.lw] = deal('none',1);

% spc.color = LD_plots.blue5;
% LD_plots.plot_pts(pts_mat_crvi(pSjh,i_crv),meta0,nse_plot1,spc);
% spc.color = LD_plots.green5;
% LD_plots.plot_pts(pts_mat_crvi(pSjhR1,i_crv),meta0,nse_plot1,spc);
% [spc.lspec,spc.lw] = deal('-',0.5);
spc.color = LD_plots.purple1;
% LD_plots.plot_pts(pts_mat_crvi(pSj0R1,i_crv),meta0,nse_plot1,spc);
% LD_plots.plot_pts(get_pts_mat(pSj0R1),meta0,nse_plot1,spc);
spc.color = LD_plots.purple5;
% LD_plots.plot_pts(pts_mat_crvi(pSj1R1,i_crv),meta0,nse_plot1,spc);
% LD_plots.plot_pts(get_pts_mat(pSj1R1),meta0,nse_plot1,spc);
spc.color = LD_plots.green4;
LD_plots.plot_pts(pts_mat_crvi(pSjR1,i_crv),meta0,nse_plot1,spc);
% LD_plots.plot_pts(get_pts_mat(pSjR1),meta0,nse_plot1,spc);

spc.color = LD_plots.green5;
LD_plots.plot_pts(pts_mat_crvi(pSjR1_1,i_crv),meta0,nse_plot1,spc);


% nse_plot1.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])

% pause
% nse_plot1.show_menubar();
nse_plot1.show_toolbar();

return

tdim_E_mesh = 6;
tdim_E_wdth = 3;
tdim_E_hght = 4;
spc.mspec = 'none';
spc.ms = 3;
spc.lspec = '-';
spc.lw = 1;
est_plot1 = LD_denoise_plots.plot_curve_estimates(nse_plot1, ...
    LD_plots('error_estimation', ...
    [tdim_E_mesh tdim_E_mesh],[tdim_E_hght tdim_E_wdth],[tdim_E_mesh 1],scrn_id), ...
    spc,Sref,Snse,i_crv);
% est_plot1.show_menubar();
est_plot1.show_toolbar();

% return

%% inspecting noisy spectrums

tdim_R_mesh = tdim_E_mesh;
tdim_R_wdth = tdim_E_wdth;
tdim_R_hght = tdim_E_hght;
rsp_plot1 = LD_denoise_plots.plot_rowspace_image(LD_plots('rowspace_image', ...
    [tdim_R_mesh tdim_R_mesh],[tdim_R_hght tdim_R_wdth],[tdim_R_mesh tdim_E_wdth+1],scrn_id), ...
    spc,Sref,Snse, ...
    Rnsvd_g,rwimg_Rn_g,Rnsvd_g_n1,rwimg_Rn_g_n1,rwimg_Rn_st_g_n1);
rsp_plot1.show_toolbar();

% rsp_plot1 = LD_plots('rowspace_image', ...
%     [tdim_R_mesh tdim_R_mesh],[tdim_R_hght tdim_R_wdth],[tdim_R_mesh tdim_E_wdth+1],1);
% % rsp_plot1.show_menubar();
% rsp_plot1.show_toolbar();
%
% plt = rsp_plot1;
% [tdim1,tdim2] = deal(3,2);
% plt = plt.init_tiles_safe(tdim1,tdim2);
% hold(plt.axs, 'on');
% box(plt.axs,'on');
% axs = plt.axs;
% axs_mat = plt.axs_mat();

% return
%{
    matrix diagnostics
%}
tdim_M_mesh = tdim_E_mesh;
tdim_M_wdth = tdim_E_wdth;
tdim_M_hght = tdim_E_hght;
mat_plot1 = LD_plots('matrix_SVD', ...
    [tdim_M_mesh tdim_M_mesh],[tdim_M_hght tdim_M_wdth],[tdim_M_mesh tdim_E_wdth+1],1);
% mat_plot1.show_menubar();
mat_plot1.show_toolbar();

plt = mat_plot1;
[tdim1,tdim2] = deal(3,2);
plt = plt.init_tiles_safe(tdim1,tdim2);
hold(plt.axs, 'on');
box(plt.axs,'on');
axs = plt.axs;
axs_mat = plt.axs_mat();

ASVD = @(c_) c_{1}*(diag(c_{2}))*(c_{3}');
r2AB = @(A_,B_) sum((A_*B_).^2,1);
rAB = @(A_,B_) sqrt(r2AB(A_,B_));

cevl = { @(c_,r_) c_{2}./max(c_{2}), @(c_,r_) 1:length(c_{2}), '\sigma_i / \sigma_1'; ...
         @(c_,r_) rAB(ASVD(r_),c_{3}), @(c_,r_) 1:length(c_{2}), '\| A \tilde{\mathbf{v}}_i \|_2'; ...
         @(c_,r_) rAB(ASVD(c_),r_{3}), @(c_,r_) (1:length(c_{2})), '\| \tilde{A} \mathbf{v}_i \|^2'  ; ...
         };
         % @(c_,r_) c_{2}, @(c_,r_) 1:length(c_{2}), '\sigma_i'; ...
         % @(c_,r_) c_{2}./max(c_{2}), @(c_,r_) 1:length(c_{2}), '\sigma_i / \sigma_1'; ...
         % @(c_,r_) sqrt(sum((r_{1}*diag(r_{2})*(r_{3}')*c_{3}).^2,1)), @(c_,r_) (1:length(c_{2})), '\| A \tilde{\mathbf{v}}_i \|_2' };
         % @(c_,r_) ( (c_{2}-r_{2}) ),   @(c_,r_) (1:length(c_{2})), '\tilde{\sigma}_i - \sigma_i' };
         % @(c_,r_) ( c_{2}(2:end)./c_{2}(1:(end-1)) ),   @(c_,r_) 1:(length(c_{2})-1), '\sigma_{i+1} / \sigma_i' };
exp_names = {'R^{(1)} [\mathrm{multinomial}]' ; 'R^{(1)} [\mathrm{Legendre}]'};
expcell{:,:,1} = {  R1svd_gf1jet.U, R1svd_gf1jet.s, R1svd_gf1jet.V ; ...
                    R1svd_gf1jet_n1.U, R1svd_gf1jet_n1.s, R1svd_gf1jet_n1.V ; ...
                };
                % Lsvd_gf1jet_n1.U, Lsvd_gf1jet_n1.s, Lsvd_gf1jet_n1.V ; ...
                % Lsvd_gf1jet_n1h.U, Lsvd_gf1jet_n1h.s, Lsvd_gf1jet_n1h.V ; ...
                % Lsvd_gf1jet_n2.U, Lsvd_gf1jet_n2.s, Lsvd_gf1jet_n2.V ; ...
                % Lsvd_gf1jet_n2h.U, Lsvd_gf1jet_n2h.s, Lsvd_gf1jet_n2h.V ; ...
                    % Lsvd_gf1jet_n2.U, Lsvd_gf1jet_n2.s, Lsvd_gf1jet_n2.V };
expcell{:,:,2} = {  Rnsvd_g.U, Rnsvd_g.s, Rnsvd_g.V ; ...
                    Rnsvd_g_n1.U, Rnsvd_g_n1.s, Rnsvd_g_n1.V ; ...
                };
                % R1svd_gf1jet.U, R1svd_gf1jet.s, R1svd_gf1jet.V ; ...
                % R1svd_gf1jet_n1.U, R1svd_gf1jet_n1.s, R1svd_gf1jet_n1.V ; ...
                % R1svd_gf1jet_n1h.U, R1svd_gf1jet_n1h.s, R1svd_gf1jet_n1h.V ; ...
                % R1svd_gf1jet_n2.U, R1svd_gf1jet_n2.s, R1svd_gf1jet_n2.V ; ...
                % R1svd_gf1jet_n2h.U, R1svd_gf1jet_n2h.s, R1svd_gf1jet_n2h.V ; ...
                    % R1svd_gf1jet_n2.U, R1svd_gf1jet_n2.s, R1svd_gf1jet_n2.V };
lcell = {'none'; 'none'; 'none'};
mcell = { 'o' ; 's'; 'v' ; '^'; '+' ; 'x'};
% cmat = [ 0 0 0 ; LD_plots.blue5 ; LD_plots.red5 ; LD_plots.orange1 ; LD_plots.green4 ; LD_plots.green5];
cmat = [ 0 0 0 ; LD_plots.red5 ; LD_plots.blue5 ; LD_plots.orange1 ; LD_plots.green4 ; LD_plots.green5];
[nevl,nexp,nnse] = deal(size(cevl,1), size(expcell,3), size(expcell{:,:,1},1));
for iexp = 1:nexp
    title(axs_mat(1,iexp), ['$$' exp_names{iexp} '$$'],'Interpreter','Latex','FontSize',10);
end
for ievl = 1:nevl
    for iexp = 1:nexp
        for inse = 1:nnse
            plot(axs_mat(ievl,iexp), ...
                cevl{ievl,2}( expcell{:,:,iexp}(inse,:), expcell{:,:,iexp}(1,:) ), ...
                cevl{ievl,1}( expcell{:,:,iexp}(inse,:), expcell{:,:,iexp}(1,:) ), ...
                'LineStyle',lcell{ievl}, ...
                'Marker',mcell{inse}, ...
                'Color',cmat(inse,:), ...
                'MarkerSize',4,'LineWidth',1)
        end
    end
    ylabel(axs_mat(ievl,:),['$$' cevl{ievl,3} '$$'],'Interpreter','Latex','FontSize',16)
end

set(axs(:), 'TickLabelInterpreter','Latex','FontSize',12);
set(axs_mat(1,:),'YScale','log', 'XScale','linear');
% set(axs_mat(1:2,:),'YScale','log', 'XScale','linear');
% set(axs_mat(3,:),'YScale','linear', 'XScale','linear');
set(axs_mat(2:3,:),'YScale','log', 'XScale','linear');
% ylabel(axs_mat(ievl,:),['$$' cevl{ievl,3} '$$'],'Interpreter','Latex','FontSize',16)
