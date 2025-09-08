clear;
close all;

scrn_id = 1;

% dir_name = '../data_directory';
% dir_name = '../denoise_data_directory';
% dir_name = '../denoise_data_directory/Gaussian_IC_perturbation';
% dir_name = '../denoise_data_directory/Gaussian_IC_perturbation/rendering_data';
dir_name = '../denoise_data_directory/Uniform_IC_perturbation/rendering_data';
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

noise_level = 0;
% noise_level = 1;
% noise_level = 2;
% Snse = Sref;
Snse = LD_observations_set(dir_name,eqn_name,['noise' num2str(noise_level)],'DoP853', dat_suff);
% solspc_nse_plot = LD_plots.plot_solspc(Snse,LD_plots('Snse',[4 4],[1 4],[2 1],1),spc);

Sref.crvs = Sref.make_curve_array;
Snse.crvs = Snse.make_curve_array;

pts_ref_cell = Sref.pts_cell();
pts_nse_cell = Snse.pts_cell();

Rsvd_g_names = { 'Rsvd_g' ; 'Rsvd_h_g' };
jet_sol_names = { '.jsol_h'; ...
'.jsol_h_Rk'; ...
'.jsol_0_Rk'; ...
'.jsol_1_Rk'; ...
};
[Rsvd_g_0,Rsvd_h_g_0,theta_jsh,pSj] = Snse.read_jet_sol_h_data(Rsvd_g_names,'.theta_mat', jet_sol_names);
pSj_cells = LD_observations_set.pts_struct_2_cell(pSj);
pSj_h_cell = pSj_cells(:,1);
pSj_h_Rk_cell = pSj_cells(:,2);
pSj_0_Rk_cell = pSj_cells(:,3);
pSj_1_Rk_cell = pSj_cells(:,4);
pSj_Rk_1_cell = LD_observations_set.combine_pts_cells( pSj_cells(:,3:4),pts_nse_cell );
dxuk_Rk = Snse.read_dxuk_data('.dxuk_Rk');
% pSj_Rk_cell = Snse.

[Rsvd_g_0_tru,Rsvd_h_g_0_tru,theta_jsh_tru,pSj_tru] = Sref.read_jet_sol_h_data(Rsvd_g_names,'.theta_mat', jet_sol_names);
pSj_cells_tru = LD_observations_set.pts_struct_2_cell(pSj_tru);
pSj_h_cell_tru = pSj_cells_tru(:,1);
pSj_h_Rk_cell_tru = pSj_cells_tru(:,2);
pSj_0_Rk_cell_tru = pSj_cells_tru(:,3);
pSj_1_Rk_cell_tru = pSj_cells_tru(:,4);
pSj_Rk_1_cell_tru = LD_observations_set.combine_pts_cells( pSj_cells_tru(:,3:4),pts_ref_cell );
dxuk_Rk_tru = Sref.read_dxuk_data('.dxuk_Rk');


% Snse1 = LD_observations_set(dir_name,eqn_name,['noise' num2str(1)],'DoP853', dat_suff);

% Snse_dns1 = LD_observations_set(dir_name,eqn_name,['noise' num2str(noise_level)],'DoP853','.jsol_Rk_1', dat_suff);

% Lsvd_gf1jet = Sref.read_LD_svd('Lsvd_global_f1jet');
% Rksvd_gf1jet = Sref.read_LD_svd('Rksvd_global_f1jet');
%     Lsvd_gf1jet_h = Sref.read_LD_svd('Lsvd_global_f1jet_h');
%     Rksvd_gf1jet_h = Sref.read_LD_svd('Rksvd_global_f1jet_h');
% Rnsvd_g = Sref.read_LD_svd('Rnsvd_global',fam_name,bor);
% rwimg_Rn_g = Sref.read_rowspace_image('Rnsvd_global',fam_name,bor);

% Lsvd_gf1jet_n1 = Snse1.read_LD_svd('Lsvd_global_f1jet');
% Rksvd_gf1jet_n1 = Snse1.read_LD_svd('Rksvd_global_f1jet');
%     Lsvd_gf1jet_n1h = Snse1.read_LD_svd('Lsvd_global_f1jet_h');
%     Rksvd_gf1jet_n1h = Snse1.read_LD_svd('Rksvd_global_f1jet_h');
% Rnsvd_g_n1 = Snse1.read_LD_svd('Rnsvd_global',fam_name,bor);
% rwimg_Rn_g_n1 = Snse1.read_rowspace_image('Rnsvd_global',fam_name,bor);
% rwimg_Rn_st_g_n1 = Snse1.read_rowspace_image('Rnsvd_strue_global',fam_name,bor);

%% inspecting curves
ncrv_ref = Sref.ncrv;

i_crv = 1;
% i_crv = 2;
% i_crv = 3;

tdim_S_mesh = 13;
tdim_S_wdth = tdim_S_mesh;
tdim_S_hght = 3;

nplt = 2;
slnspc_plts = LD_plots.empty(nplt,0);
for i = 1:nplt
    slnspc_plts(i) = LD_plots(['Sdns_plt', num2str(i)], ...
                    [tdim_S_mesh tdim_S_mesh],...
                    [tdim_S_hght tdim_S_wdth],[(tdim_S_hght+1)*(i+1) 1], ...
                scrn_id);
end
slnspc_ref_plt = slnspc_plts(1);
spc.color = [0 0 0 0.25];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'.',5);
slnspc_ref_plt = LD_plots.plot_pts(pts_ref_cell(1:ncrv_ref), ...
                            meta0,slnspc_ref_plt,spc);
spc.color = [0 0 0 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,'o',2);
slnspc_ref_plt = LD_plots.plot_pts(pts_ref_cell(i_crv), ...
                        meta0,slnspc_ref_plt,spc);
spc.color = [1 0 0 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,'s',2);
slnspc_ref_plt = LD_plots.plot_pts(pts_nse_cell(i_crv), ...
                        meta0,slnspc_ref_plt,spc);
spc.color = [LD_plots.green1 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'d',2);
slnspc_ref_plt = LD_plots.plot_pts(pSj_h_cell(i_crv), ...
                        meta0,slnspc_ref_plt,spc);
spc.color = [LD_plots.green4 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'d',2);
slnspc_ref_plt = LD_plots.plot_pts(pSj_h_Rk_cell(i_crv), ...
                        meta0,slnspc_ref_plt,spc);
spc.color = [LD_plots.purple1 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'<',2);
slnspc_ref_plt = LD_plots.plot_pts(pSj_0_Rk_cell(i_crv), ...
                    meta0,slnspc_ref_plt,spc);
spc.color = [LD_plots.purple5 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'>',2);
slnspc_ref_plt = LD_plots.plot_pts(pSj_1_Rk_cell(i_crv), ...
                    meta0,slnspc_ref_plt,spc);
spc.color = [0 0 1 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,'s',2);
slnspc_ref_plt = LD_plots.plot_pts(pSj_Rk_1_cell(i_crv), ...
                        meta0,slnspc_ref_plt,spc);
% slnspc_ref_plt.show_menubar();
slnspc_ref_plt.show_toolbar();


slnspc_nse_plt = slnspc_plts(2);
spc.color = [1 0 0 0.25];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'*',2);
slnspc_nse_plt = LD_plots.plot_pts(pts_nse_cell(1:ncrv_ref), ...
                            meta0,slnspc_nse_plt,spc);
spc.color = [1 0 0 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,'s',4);
slnspc_nse_plt = LD_plots.plot_pts(pts_nse_cell(i_crv), ...
                        meta0,slnspc_nse_plt,spc);
spc.color = [0 0 0 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,'o',2);
slnspc_nse_plt = LD_plots.plot_pts(pts_ref_cell(i_crv), ...
                        meta0,slnspc_nse_plt,spc);
spc.color = [LD_plots.green1 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'d',2);
slnspc_nse_plt = LD_plots.plot_pts(pSj_h_cell(i_crv), ...
                        meta0,slnspc_nse_plt,spc);
spc.color = [LD_plots.green4 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'d',2);
slnspc_nse_plt = LD_plots.plot_pts(pSj_h_Rk_cell(i_crv), ...
                        meta0,slnspc_nse_plt,spc);
spc.color = [LD_plots.purple1 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'<',2);
slnspc_nse_plt = LD_plots.plot_pts(pSj_0_Rk_cell(i_crv), ...
                    meta0,slnspc_nse_plt,spc);
spc.color = [LD_plots.purple5 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('none',1,'>',2);
slnspc_nse_plt = LD_plots.plot_pts(pSj_1_Rk_cell(i_crv), ...
                    meta0,slnspc_nse_plt,spc);
spc.color = [0 0 1 1];
[spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,'s',4);
slnspc_nse_plt = LD_plots.plot_pts(pSj_Rk_1_cell(i_crv), ...
                        meta0,slnspc_nse_plt,spc);
% slnspc_nse_plt.show_menubar();
slnspc_nse_plt.show_toolbar();




i_crv = 1;
spc.color = [0 0 0 0.2];
plt0 = LD_plots('Sdns', ...
                [tdim_S_mesh tdim_S_mesh],...
                    [tdim_S_hght tdim_S_wdth],[tdim_S_hght 1], ...
                scrn_id);
nse_plot1 = LD_denoise_plots.plot_denoised_trajectories(plt0, ...
                                                        spc, ...
                                                        Sref,Snse,i_crv);
% nse_plot1.show_menubar();
nse_plot1.show_toolbar();
% return


return













[Rsvd_g_0,Rsvd_h_g_0,theta_jsh,pSjh,pSjhRk,pSj0Rk,pSj1Rk] = Snse.read_jet_sol_h_data(Rsvd_g_names,'.theta_mat', jet_sol_names);
% [theta_jsh,pSjh,pSjhRk,pSj0Rk,pSj1Rk,pSjRk] = Snse.read_jet_sol_h_data('.theta_mat', jet_sol_names);
% [theta_jsh_1,pSjh_1,pSjhRk_1,pSj0Rk_1,pSj1Rk_1,pSjRk_1] = Snse_dns1.read_jet_sol_h_data('.theta_mat', jet_sol_names);

get_pts_mat = @(pS_) reshape(pS_.pts_in,meta0.ndim,[]);
pts_mat_crvi = @(pS_,i_) reshape(pS_.pts_in( (pS_.pts_crv_inds(1,i_)):(pS_.pts_crv_inds(2,i_)) ) , meta0.ndim,[]);

[spc.mspec,spc.ms] = deal('.',5);
[spc.lspec,spc.lw] = deal('none',1);

% spc.color = LD_plots.blue5;
% LD_plots.plot_pts(pts_mat_crvi(pSjh,i_crv),meta0,nse_plot1,spc);
% spc.color = LD_plots.green5;
% LD_plots.plot_pts(pts_mat_crvi(pSjhRk,i_crv),meta0,nse_plot1,spc);
% [spc.lspec,spc.lw] = deal('-',0.5);
% spc.color = LD_plots.purple1;
% LD_plots.plot_pts(pts_mat_crvi(pSj0Rk,i_crv),meta0,nse_plot1,spc);
% % LD_plots.plot_pts(get_pts_mat(pSj0Rk),meta0,nse_plot1,spc);
% spc.color = LD_plots.purple5;
% LD_plots.plot_pts(pts_mat_crvi(pSj1Rk,i_crv),meta0,nse_plot1,spc);
% % LD_plots.plot_pts(get_pts_mat(pSj1Rk),meta0,nse_plot1,spc);

% spc.color = LD_plots.green5;
% LD_plots.plot_pts(pts_mat_crvi(pSjRk_1,i_crv),meta0,nse_plot1,spc);


% nse_plot1.write_figure('png',[getenv('HOME') '/Desktop/MATLAB_OUTPUT/'])

% pause

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
expcell{:,:,1} = {  Rksvd_gf1jet.U, Rksvd_gf1jet.s, Rksvd_gf1jet.V ; ...
                    Rksvd_gf1jet_n1.U, Rksvd_gf1jet_n1.s, Rksvd_gf1jet_n1.V ; ...
                };
                % Lsvd_gf1jet_n1.U, Lsvd_gf1jet_n1.s, Lsvd_gf1jet_n1.V ; ...
                % Lsvd_gf1jet_n1h.U, Lsvd_gf1jet_n1h.s, Lsvd_gf1jet_n1h.V ; ...
                % Lsvd_gf1jet_n2.U, Lsvd_gf1jet_n2.s, Lsvd_gf1jet_n2.V ; ...
                % Lsvd_gf1jet_n2h.U, Lsvd_gf1jet_n2h.s, Lsvd_gf1jet_n2h.V ; ...
                    % Lsvd_gf1jet_n2.U, Lsvd_gf1jet_n2.s, Lsvd_gf1jet_n2.V };
expcell{:,:,2} = {  Rnsvd_g.U, Rnsvd_g.s, Rnsvd_g.V ; ...
                    Rnsvd_g_n1.U, Rnsvd_g_n1.s, Rnsvd_g_n1.V ; ...
                };
                % Rksvd_gf1jet.U, Rksvd_gf1jet.s, Rksvd_gf1jet.V ; ...
                % Rksvd_gf1jet_n1.U, Rksvd_gf1jet_n1.s, Rksvd_gf1jet_n1.V ; ...
                % Rksvd_gf1jet_n1h.U, Rksvd_gf1jet_n1h.s, Rksvd_gf1jet_n1h.V ; ...
                % Rksvd_gf1jet_n2.U, Rksvd_gf1jet_n2.s, Rksvd_gf1jet_n2.V ; ...
                % Rksvd_gf1jet_n2h.U, Rksvd_gf1jet_n2h.s, Rksvd_gf1jet_n2h.V ; ...
                    % Rksvd_gf1jet_n2.U, Rksvd_gf1jet_n2.s, Rksvd_gf1jet_n2.V };
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
