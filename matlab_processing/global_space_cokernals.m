clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

xrange = 0;
% xrange = 1;
% xrange = 2;

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

bor = 10;
% bor = 9;
% bor = 8;
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
solspc_ref_plot = LD_plots.plot_solspc(Sref,LD_plots('Sref',[5 5],[1 5],[1 1],1),spc);

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

k_clst = 2:5;

%%%%%%%%%%%%%%%%%%%%%%%%

% Ytns_L_check = AYLrstsvd.Ytns_L
% Ytns_L_check = Sref_Lsvd.Vtns(:,1:(min(Sref_Lsvd.rvec)),:);
Ytns_L_check = Sref_Lsvd.Vtns(:,1:(max(Sref_Lsvd.rvec)),:);
cmat_Y_L = LD_aux.compute_frobenius_closeness_matrix(Ytns_L_check);
dmat_Y_L = 1.0 - cmat_Y_L;

SC_YL_vec = nan(size(k_clst));
for ik = 1:length(k_clst)
    SC_YL_vec(ik) = getfield(LD_aux.compute_medoid_silhouette(LD_aux.naive_k_medoid(dmat_Y_L,k_clst(ik))), ...
                                'silhouette_coeff');
end
k_clst_YL_silhouettes = [k_clst; SC_YL_vec]
[~,i_SC_YL_max] = max(SC_YL_vec);

k_clusters_YL = k_clst(i_SC_YL_max);
[medoid_YL_pckg,greedy_medoid_YL_pckg,improved_on_greedy_YL] = LD_aux.naive_k_medoid(dmat_Y_L,k_clusters_YL)
[~,medoid_YL_pckg_detailed] = LD_aux.post_process_medoid_package(medoid_YL_pckg)


% Ktns_A = Asvd.Vtns(:,(max(Asvd.rvec)+1):end,:);
Ktns_A = Asvd.Vtns(:,(min(Asvd.rvec)+1):end,:);
cmat_K_A = LD_aux.compute_frobenius_closeness_matrix(Ktns_A);
dmat_K_A = 1.0 - cmat_K_A;

SC_KA_vec = nan(size(k_clst));
for ik = 1:length(k_clst)
    SC_KA_vec(ik) = getfield(LD_aux.compute_medoid_silhouette(LD_aux.naive_k_medoid(dmat_K_A,k_clst(ik))), ...
                                'silhouette_coeff');
end
k_clst_KA_silhouettes = [k_clst; SC_KA_vec]
[~,i_SC_KA_max] = max(SC_KA_vec);

k_clusters_KA = k_clst(i_SC_KA_max);
[medoid_KA_pckg,greedy_medoid_KA_pckg,improved_on_greedy_KA] = LD_aux.naive_k_medoid(dmat_K_A,k_clusters_KA)
[~,medoid_KA_pckg_detailed] = LD_aux.post_process_medoid_package(medoid_KA_pckg)


Ktns_A_YLrst = LD_aux.compute_Ktns_A_YLrst(AYLrstsvd);
cmat_K_A_YLrst = LD_aux.compute_frobenius_closeness_matrix(Ktns_A_YLrst);
dmat_K_A_YLrst = 1.0 - cmat_K_A_YLrst;

SC_KAYL_vec = nan(size(k_clst));
for ik = 1:length(k_clst)
    SC_KAYL_vec(ik) = getfield(LD_aux.compute_medoid_silhouette(LD_aux.naive_k_medoid(dmat_K_A_YLrst,k_clst(ik))), ...
                                'silhouette_coeff');
end
k_clst_KAYL_silhouettes = [k_clst; SC_KAYL_vec]
[~,i_SC_KAYL_max] = max(SC_KAYL_vec);

k_clusters_KAYL = k_clst(i_SC_KAYL_max);
[medoid_KAYL_pckg,greedy_medoid_KAYL_pckg,improved_on_greedy_KAYL] = LD_aux.naive_k_medoid(dmat_K_A_YLrst,k_clusters_KAYL)
[~,medoid_KAYL_pckg_detailed] = LD_aux.post_process_medoid_package(medoid_KAYL_pckg)

%%%%%%%%%%%%%%%%%%%%%%%

% [medoid_pckg_check,k_clst_silh_check] = deal(medoid_YL_pckg,k_clst_YL_silhouettes);
[medoid_pckg_check,k_clst_silh_check] = deal(medoid_KA_pckg,k_clst_KA_silhouettes);
% [medoid_pckg_check,k_clst_silh_check] = deal(medoid_KAYL_pckg,k_clst_KAYL_silhouettes);

%%%%%%%%%%%%%%%%%%%%%%%

[icrv_med,ncrv_clst,net_dfro,clst_membership,local_dfro,med2med_dfro,dmat_check] = LD_aux.unpack_medoid_package(medoid_pckg_check);
clst_info = LD_aux.get_cluster_info(medoid_pckg_check);
k_clst_silh_check

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
    solspc_ref_plot_i = LD_plots(['Sref_clst' num2str(i)],[5 5],[1 5],[1+i 1],1);

    spc.color = col_med_plot_mat(i,:);
    spc.lw = 2;
    solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_med_plot(i));
    solspc_ref_plot_i = LD_plots.plot_solspc(Sref,solspc_ref_plot_i,spc,icrv_med_plot(i));

    spc.color = [col_med_plot_mat(i,:) alpha_cluster];
    spc.lw = 1;
    solspc_ref_plot_i = LD_plots.plot_solspc(Sref,solspc_ref_plot_i,spc,icrv_vec(clst_membership == isrt_ncrv_clst(i)));
end
solspc_plot_init_lims = solspc_ref_plot.get_axis_lims;

[xlim_hist,ylim_hist] = deal(nan(length(nmed_plot),1));
aux_fig1 = LD_plots('aux1',[6 6],[1 6],[6 1],1);
aux_fig1 = aux_fig1.init_tiles_safe(1,nmed_plot);
hold(aux_fig1.axs, 'on');
box(aux_fig1.axs,'on');
axs_aux1 = aux_fig1.axs;
for i = 1:nmed_plot
    histogram(axs_aux1(i), dmat_check(:,icrv_med_plot(i)), 'FaceColor', col_med_plot_mat(i,:));
    xlim_hist(i) = max(get(axs_aux1(i),'XLim'));
    ylim_hist(i) = max(get(axs_aux1(i),'YLim'));
end
set(axs_aux1,'XLim',[0,max(dmat_check(:))],'XLim',[0,max(xlim_hist)],'YLim',[0,max(ylim_hist)]);
