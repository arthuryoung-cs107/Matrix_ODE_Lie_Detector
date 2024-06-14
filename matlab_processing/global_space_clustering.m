clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

xrange = 0;
% xrange = 1;
% ode_name = 'Duffing';
ode_name = 'VanDerPol';

eqn_name = [ode_name '_xrange' num2str(xrange)];
% eqn_name = [ode_name '_xrange' num2str(xrange) '_extrap'];

fam_name = 'Chebyshev1';
% fam_name = 'Legendre';
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
Sref_OGsvd = Sref.read_mat_svd_package(fam_name,bor,'OGmat');

Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Psvd, ...
                            Sref_Qsvd, ...
                            Sref_Gsvd, ...
                            Sref_OGsvd];

Sref_svd_glb_basic_array = LD_aux.make_global_svd_package(Sref_svd_basic_array);
Sref_Lsvd_glb = Sref_svd_glb_basic_array(1);
Sref_Rsvd_glb = Sref_svd_glb_basic_array(2);
Sref_Osvd_glb = Sref_svd_glb_basic_array(3);
Sref_Psvd_glb = Sref_svd_glb_basic_array(4);
Sref_Qsvd_glb = Sref_svd_glb_basic_array(5);
Sref_Gsvd_glb = Sref_svd_glb_basic_array(6);
Sref_OGsvd_glb = Sref_svd_glb_basic_array(7);

% [Sref_OYLglbsvd,Sref_OYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_Osvd,Sref_Lsvd_glb,1);
% [Sref_GYLglbsvd,Sref_GYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_Gsvd,Sref_Lsvd_glb,1);
% [Sref_OGYLglbsvd,Sref_OGYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_OGsvd,Sref_Lsvd_glb,1);

% Sref_svd_Lglb_reg_array =   [   Sref_OYLglbsvd, Sref_GYLglbsvd, Sref_OGYLglbsvd ];
% Sref_svd_glb_reg_array =    [   Sref_OYLglbsvd_glb, Sref_GYLglbsvd_glb, Sref_OGYLglbsvd_glb];
Sref_svd_Lglb_reg_array =   [];
Sref_svd_glb_reg_array =    [];

Sref_svd_glb_array = [Sref_svd_glb_basic_array, Sref_svd_glb_reg_array];

Sref_glb_svd_plot = LD_plots.plot_global_svds(Sref_svd_glb_array,LD_plots('svd',[7 7],[2 7],[5 1],1));
Sref_glb_svd_plot.show_toolbar();

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

Sref_svd_array = [Sref_svd_basic_array, Sref_svd_Lglb_reg_array, Sref_svd_YL_in_array, Sref_svd_YL_array];
Sref_crv_svd_plot = LD_plots.plot_curve_svds(Sref_svd_array,LD_plots('svd',[7 7],[2 7],[7 1],1));

%%%%%%%%%%%%%%%%%%%%%%%%

Lsvd_glb = Sref_Lsvd_glb;
Asvd = Sref_Osvd;
AYLsvd = Sref_OYLsvd;
AYLrstsvd = Sref_OYLsvd_rst;
% AYLglbsvd = Sref_OYLglbsvd_glb;

% Lsvd_glb = Sref_Lsvd_glb;
% Asvd = Sref_Gsvd;
% AYLsvd = Sref_GYLsvd;
% AYLrstsvd = Sref_GYLsvd_rst;
% % AYLglbsvd = Sref_GYLglbsvd_glb;

% Lsvd_glb = Sref_Lsvd_glb;
% Asvd = Sref_OGsvd;
% AYLsvd = Sref_OGYLsvd;
% AYLrstsvd = Sref_OGYLsvd_rst;
% % AYLglbsvd = Sref_OGYLglbsvd_glb;


[ndep,nvar] = deal(meta0.ndep,meta0.ndep+1);
[ncrv,nrow,ncol] = deal(length(Asvd.rvec),Asvd.nrow,Asvd.ncol);
% ncol_glb = AYLglbsvd.ncol;
% kappa_L = (ncol-ncol_glb)/nvar;
% rho_L = Lsvd_glb.ncol - kappa_L;

Amat = Asvd.matT';
Atns = permute(reshape(Amat',ncol,[],ncrv),[2 1 3]);
Ytns_L = AYLrstsvd.Ytns_L;
AtnsYL_rst = LD_aux.Atns_Ytns_mult(Atns,Ytns_L);

%%%%%%%%%%%%%%%%%%%%%%%%

svd_datum = Sref_Gsvd;
[~,isort_rvec_datum] = sort(svd_datum.rvec);

icrv_centre = 1;
% icrv_centre = isort_rvec_datum(1) % lowest rank
% icrv_centre = isort_rvec_datum(ceil(ncrv/2)) % median rank
% icrv_centre = isort_rvec_datum(end) % max rank

[rho_PL,ncol_AYL,rho_AYL] = deal(AYLrstsvd.rho_PL,AYLrstsvd.ncol,max(AYLrstsvd.rvec))
kappa_AYL = ncol_AYL - rho_AYL;

Ktns_A_YLrst = LD_aux.Ytns_Vtns_mult(Ytns_L,AYLrstsvd.Vtns(:,(rho_AYL+1):end,:)); % AYLmat kernals of each curve
Kmat_A_YLrst = reshape(Ktns_A_YLrst,ncol,[]);

KTK_A_YLrst = mat2cell(Kmat_A_YLrst'*Kmat_A_YLrst,kappa_AYL*ones(1,ncrv), kappa_AYL*ones(1,ncrv));
K_A_YLrst = cell(ncrv);
for i = 1:ncrv
    for j = 1:ncrv
        K_A_YLrst{i,j} = Ktns_A_YLrst(:,:,i)*KTK_A_YLrst{i,j};
    end
end
K_A_YLrst = mat2cell(Kmat_A_YLrst'*Kmat_A_YLrst,kappa_AYL*ones(1,ncrv), kappa_AYL*ones(1,ncrv));

mflag_lw = logical(tril(ones(ncrv),-1)); % flags of lower triangular elements
vflag_lw = reshape(mflag_lw,[],1); % column stack of lower triangular flags
vinds_full = (1:(ncrv*ncrv))'; % 1 thru ncrv^2
vinds_full_mat = reshape(vinds_full,ncrv,ncrv); % matrix of column stacked indices
vinds_lw = vinds_full(vflag_lw); % column stacked indices of strictly lower triangular elements

KTK_A_YLrst_lwvec = KTK_A_YLrst(vflag_lw); % column stacked sub cells for lower triangular elements

% fro_nrmmag = @(mat_) norm(mat_,'fro')/sqrt(size(mat_,2));
fro_nrmmag = @(mat_) norm(mat_,'fro')/sqrt(size(mat_,2)-1);
% avg_colmag = @(mat_) mean(sqrt(sum(mat_.*mat_,1)),2);
avg_colmag = @(mat_) sum( sqrt(sum(mat_.*mat_,1)), 2 )/( size(mat_,2)-1 );
avg_rowmag = @(mat_) avg_colmag(mat_');

avg_fmag_KTK_A_YLrst_lwvec = cellfun(fro_nrmmag,KTK_A_YLrst_lwvec);
avg_cmag_KTK_A_YLrst_lwvec = cellfun(avg_colmag,KTK_A_YLrst_lwvec);
avg_rmag_KTK_A_YLrst_lwvec = cellfun(avg_rowmag,KTK_A_YLrst_lwvec);

evl_mag_KTK_A_YLrst_lwvec = avg_fmag_KTK_A_YLrst_lwvec;
% evl_mag_KTK_A_YLrst_lwvec = avg_cmag_KTK_A_YLrst_lwvec;
% evl_mag_KTK_A_YLrst_lwvec = avg_rmag_KTK_A_YLrst_lwvec;

[~,isrt_evl_mag_KTK_A_YLrst_lwvec] = sort(evl_mag_KTK_A_YLrst_lwvec); % sorted point stat of cokernal magnitudes
[idst_evl_lwvec,ingh_evl_lwvec] = deal(isrt_evl_mag_KTK_A_YLrst_lwvec(1),isrt_evl_mag_KTK_A_YLrst_lwvec(end))
[idst_evl,ingh_evl] = deal(vinds_lw(idst_evl_lwvec),vinds_lw(ingh_evl_lwvec))
[icrv_dst1,icrv_dst2] = find(vinds_full_mat==idst_evl)
[icrv_ngh1,icrv_ngh2] = find(vinds_full_mat==ingh_evl)
evl_mag_distant_lwvec = evl_mag_KTK_A_YLrst_lwvec(idst_evl_lwvec)
evl_mag_neighbor_lwvec = evl_mag_KTK_A_YLrst_lwvec(ingh_evl_lwvec)

KTK_A_YLrst_evl_mag_lwmat = zeros(ncrv);
KTK_A_YLrst_evl_mag_lwmat(vflag_lw) = evl_mag_KTK_A_YLrst_lwvec;
KTK_A_YLrst_evl_mag = KTK_A_YLrst_evl_mag_lwmat + eye(ncrv) + KTK_A_YLrst_evl_mag_lwmat';

[min_evlmags,med_evlmags,avg_evlmags,max_evlmags] = LD_aux.mat_colstats(KTK_A_YLrst_evl_mag);

evlstat_evlmags = avg_evlmags;

[~,isrt_evlstat_evlmags] = sort(evlstat_evlmags);
[icrv_loneliest,icrv_frndliest] = deal(isrt_evlstat_evlmags(1),isrt_evlstat_evlmags(end))
% [icrv_loneliest,icrv_frndliest] = deal(isrt_evlstat_evlmags(2),isrt_evlstat_evlmags(end))

k_clusters = 2;
[icrv_med,ncrv_clst,net_dfro,clst_membership,local_dfro,med2med_dfro] = LD_aux.naive_k_medoid(1.0-KTK_A_YLrst_evl_mag,k_clusters)

%%%%%%%%%%%%%%%%%%%%%%%%

spc.lw = 2.0;

% spc.color = LD_plots.orange1;
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_dst1);
% spc.color = LD_plots.orange2;
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_dst2);
%
% spc.color = LD_plots.blue1;
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_ngh1);
% spc.color = LD_plots.blue2;
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_ngh2);
%
spc.color = LD_plots.purple1;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_loneliest);
spc.color = LD_plots.purple5;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_frndliest);

col_med_plot_mat = LD_plots.greens;

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

for i = 1:nmed_plot
    spc.color = col_med_plot_mat(i,:);
    solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_med_plot(i));
    % pause
end

solspc_plot_init_lims = solspc_ref_plot.get_axis_lims;

aux_fig1 = LD_plots('aux1',[5 5],[1 5],[5 1],1);
aux_fig1 = aux_fig1.init_tiles_safe(1,4);
hold(aux_fig1.axs, 'on');
box(aux_fig1.axs,'on');
axs_aux1 = aux_fig1.axs;

[~,isrt_pltmg] = sort(evl_mag_KTK_A_YLrst_lwvec);

plot(axs_aux1(1),1:size(avg_cmag_KTK_A_YLrst_lwvec,1),avg_cmag_KTK_A_YLrst_lwvec(isrt_pltmg), ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 0 1]);
plot(axs_aux1(1),1:size(avg_rmag_KTK_A_YLrst_lwvec,1),avg_rmag_KTK_A_YLrst_lwvec(isrt_pltmg), ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[1 0 0]);
plot(axs_aux1(1),1:size(avg_fmag_KTK_A_YLrst_lwvec,1),avg_fmag_KTK_A_YLrst_lwvec(isrt_pltmg), ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 1 0]);

plot(axs_aux1(2),avg_fmag_KTK_A_YLrst_lwvec,avg_cmag_KTK_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 0 1]);
plot(axs_aux1(2),avg_fmag_KTK_A_YLrst_lwvec,avg_rmag_KTK_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[1 0 0]);

plot(axs_aux1(3),avg_rmag_KTK_A_YLrst_lwvec,avg_cmag_KTK_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 0 1]);
plot(axs_aux1(3),avg_rmag_KTK_A_YLrst_lwvec,avg_fmag_KTK_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 1 0]);

plot(axs_aux1(4),avg_cmag_KTK_A_YLrst_lwvec,avg_rmag_KTK_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[1 0 0]);
plot(axs_aux1(4),avg_cmag_KTK_A_YLrst_lwvec,avg_fmag_KTK_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 1 0]);

aux_fig1.show_toolbar();

return
