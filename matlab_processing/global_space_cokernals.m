clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

xrange = 0;
% xrange = 1;

ode_name = 'Duffing';
% ode_name = 'VanDerPol';
% ode_name = 'Pendulum';
% ode_name = 'Bessel';
% ode_name = 'Riccati';
% ode_name = 'Brusselator';

eqn_name = [ode_name '_xrange' num2str(xrange)];
% eqn_name = [ode_name '_xrange' num2str(xrange) '_extrap'];

fam_name = 'Chebyshev1';
% fam_name = 'Legendre';

% bor = 10;
% bor = 9;
bor = 8;
% bor = 7;
% bor = 6;
% bor = 5;
% bor = 4;
% bor = 3;
% bor = 2;

spc = LD_plots.make_default_plot_specs;

Sref = LD_observations_set(dir_name,eqn_name,'true','DoP853', dat_suff);
meta0 = Sref.meta_data;

spc.color = [0 0 0 0.2];
solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[4 4],[1 4],[1 1],1),spc);

spc.color = [1 0 0];
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,68);

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

[Sref_RYLsvd_in,Sref_RYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Rmat','YL','Lmat');
Sref_svd_YL_in_array = [    Sref_RYLsvd_in];

[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
[Sref_PYLsvd,Sref_PYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Psvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
[Sref_OGYLsvd,Sref_OGYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_OGsvd,Sref_Lsvd);

Sref_svd_YL_array = [       Sref_OYLsvd, ...
                            Sref_PYLsvd, ...
                            Sref_GYLsvd, ...
                            Sref_OGYLsvd];

Sref_svd_array = [Sref_svd_basic_array, Sref_svd_YL_in_array, Sref_svd_YL_array];
Sref_crv_svd_plot = LD_plots.plot_curve_svds(Sref_svd_array,LD_plots('svd',[9 9],[2 9],[9 1],1));

%%%%%%%%%%%%%%%%%%%%%%%%

% Asvd = Sref_Osvd;
% AYLsvd = Sref_OYLsvd;
% AYLrstsvd = Sref_OYLsvd_rst;

% Asvd = Sref_Gsvd;
% AYLsvd = Sref_GYLsvd;
% AYLrstsvd = Sref_GYLsvd_rst;

Asvd = Sref_OGsvd;
AYLsvd = Sref_OGYLsvd;
AYLrstsvd = Sref_OGYLsvd_rst;


[ndep,nvar] = deal(meta0.ndep,meta0.ndep+1);
[ncrv,nrow,ncol] = deal(length(Asvd.rvec),Asvd.nrow,Asvd.ncol);

Amat = Asvd.matT';
Atns = permute(reshape(Amat',ncol,[],ncrv),[2 1 3]);
Ytns_L = AYLrstsvd.Ytns_L;
AtnsYL_rst = LD_aux.Atns_Ytns_mult(Atns,Ytns_L);

%%%%%%%%%%%%%%%%%%%%%%%%

Ktns_A_YLrst = LD_aux.compute_Ktns_A_YLrst(AYLrstsvd);
cmat_K_A_YLrst = LD_aux.compute_frobenius_closeness_matrix(Ktns_A_YLrst);
dmat_K_A_YLrst = 1.0 - cmat_K_A_YLrst;

k_clst = 2:10;

SC_KAYL_vec = nan(size(k_clst));
for ik = 1:length(k_clst)
    k_med_pckg_i = LD_aux.naive_k_medoid(dmat_K_A_YLrst,k_clst(ik));
    k_silh_i = LD_aux.compute_medoid_silhouette(k_med_pckg_i,dmat_K_A_YLrst);
    SC_KAYL_vec(ik) = k_silh_i.silhouette_coeff;
end
k_clst_KAYL_silhouettes = [k_clst; SC_KAYL_vec]

k_clusters_KAYL = 3;
[medoid_KAYL_pckg,greedy_medoid_KAYL_pckg,improved_on_greedy_KAYL] = LD_aux.naive_k_medoid(dmat_K_A_YLrst,k_clusters_KAYL)
[~,medoid_KAYL_pckg_detailed] = LD_aux.post_process_medoid_package(medoid_KAYL_pckg,dmat_K_A_YLrst)

cmat_Y_L = LD_aux.compute_frobenius_closeness_matrix(Ytns_L);
dmat_Y_L = 1.0 - cmat_Y_L;

SC_YL_vec = nan(size(k_clst));
for ik = 1:length(k_clst)
    k_med_pckg_i = LD_aux.naive_k_medoid(dmat_Y_L,k_clst(ik));
    k_silh_i = LD_aux.compute_medoid_silhouette(k_med_pckg_i,dmat_Y_L);
    SC_YL_vec(ik) = k_silh_i.silhouette_coeff;
end
k_clst_YL_silhouettes = [k_clst; SC_YL_vec]

k_clusters_YL = k_clusters_KAYL;
[medoid_YL_pckg,greedy_medoid_YL_pckg,improved_on_greedy_YL] = LD_aux.naive_k_medoid(dmat_Y_L,k_clusters_YL)
[~,medoid_YL_pckg_detailed] = LD_aux.post_process_medoid_package(medoid_YL_pckg,dmat_Y_L)

[icrv_med,ncrv_clst,net_dfro,clst_membership,local_dfro,med2med_dfro] = LD_aux.unpack_medoid_package(medoid_KAYL_pckg);
clst_info = LD_aux.get_cluster_info(medoid_KAYL_pckg);
dmat_check = dmat_K_A_YLrst;

% [icrv_med,ncrv_clst,net_dfro,clst_membership,local_dfro,med2med_dfro] = LD_aux.unpack_medoid_package(medoid_YL_pckg);
% clst_info = LD_aux.get_cluster_info(medoid_YL_pckg);
% dmat_check = dmat_Y_L;

iicrv_check1 = 1;
icrv_check1 = icrv_med(iicrv_check1);
[~,iicrv_check2] = min(clst_info{iicrv_check1}(2,:));
% [~,iicrv_check2] = max(clst_info{iicrv_check1}(2,:));
icrv_check2 = clst_info{iicrv_check1}(1,iicrv_check2);
[K_check1,K_check2] = deal(Ktns_A_YLrst(:,:,icrv_check1),Ktns_A_YLrst(:,:,icrv_check2));
% [D_K_check12,C_K_check12] = LD_aux.compute_Zassenhaus_bases(K_check1,K_check2)
% [R_K_check12,P_K_check12] = LD_aux.compute_Zassenhaus_bases(K_check1,K_check2)
% [zas_pckg1,zas_pckg2] = LD_aux.compute_Zassenhaus_bases(K_check1,K_check2)
% rank([K_check1,K_check2])
% cond([K_check1,K_check2])

% U_check = [1 -1 0 1; 0 0 1 -1]';
% W_check = [5 0 -3 3; 0 5 -3 -2]';
% % X_check = [U_check', U_check'; W_check', zeros(size(W_check))'];
% % X_check = [U_check', U_check'; U_check', U_check'; W_check', zeros(size(W_check))'];
% % X_check = [U_check', U_check'; zeros(size(U_check))', zeros(size(U_check))' ; W_check', zeros(size(W_check))'];
% X_check = [U_check', U_check'; W_check', zeros(size(W_check))'; zeros(2,2*size(U_check,1))];
% [~,U_X_check] = lu(X_check)
% % R_X_check = rref(X_check)

% return

%%%%%%%%%%%%%%%%%%%%%%%%

spc.lw = 2.0;
col_med_plot_mat = LD_plots.greens;
alpha_cluster = 0.10;

% icrv_med_plot = icrv_med;
[ncrv_clst_sort,isrt_ncrv_clst] = sort(ncrv_clst);
% [~,isrt_ncrv_clst] = sort(ncrv_clst,'descend');
icrv_med_ncrv_sort = icrv_med(isrt_ncrv_clst)
ncrv_clst_sort

icrv_med_plot = icrv_med_ncrv_sort;

if (length(icrv_med_plot) == size(col_med_plot_mat,1))
    icrv_med_plot = icrv_med_ncrv_sort;
else
    if (length(icrv_med_plot) < size(col_med_plot_mat,1))
        col_med_plot_mat = col_med_plot_mat((size(col_med_plot_mat,1) - length(icrv_med_plot) +1):end,:);
    else
        icrv_med_plot = icrv_med_ncrv_sort((end - (size(col_med_plot_mat,1)-1)):end)
    end
end

nmed_plot = min([length(icrv_med_plot), size(col_med_plot_mat,1)]);

icrv_vec = 1:ncrv;
for i = 1:nmed_plot
    % solspc_ref_plot_i = solspc_ref_plot;
    solspc_ref_plot_i = LD_plots(['Sref_clst' num2str(i)],[4 4],[1 4],[1+i 1],1);

    spc.color = col_med_plot_mat(i,:);
    spc.lw = 2;
    solspc_ref_plot_i = LD_plots.plot_solspc(Sref,solspc_ref_plot_i,spc,icrv_med_plot(i));

    spc.color = [col_med_plot_mat(i,:) alpha_cluster];
    spc.lw = 1;
    solspc_ref_plot_i = LD_plots.plot_solspc(Sref,solspc_ref_plot_i,spc,icrv_vec(clst_membership == isrt_ncrv_clst(i)));
end
solspc_plot_init_lims = solspc_ref_plot.get_axis_lims;

aux_fig1 = LD_plots('aux1',[5 5],[1 5],[1 1],1);
aux_fig1 = aux_fig1.init_tiles_safe(1,nmed_plot);
hold(aux_fig1.axs, 'on');
box(aux_fig1.axs,'on');
axs_aux1 = aux_fig1.axs;
for i = 1:nmed_plot
    % histogram(axs_aux1(i), dmat_K_A_YLrst(:,icrv_med_plot(i)), 'FaceColor', col_med_plot_mat(i,:));
    histogram(axs_aux1(i), dmat_check(:,icrv_med_plot(i)), 'FaceColor', col_med_plot_mat(i,:));
end
% set(axs_aux1,'XLim',[0,max(dmat_K_A_YLrst(:))])
set(axs_aux1,'XLim',[0,max(dmat_check(:))])
