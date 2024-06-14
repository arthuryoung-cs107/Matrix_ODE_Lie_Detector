clear;
close all;

dir_name = '../data_directory';
dat_suff = 'lddat';

xrange = 0;
% xrange = 1;
ode_name = 'Duffing';

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
% Sref_OPsvd = Sref.read_mat_svd_package(fam_name,bor,'OPmat');
% Sref_OGsvd = Sref.read_mat_svd_package(fam_name,bor,'OGmat');

Sref_svd_basic_array = [    Sref_Lsvd, ...
                            Sref_Rsvd, ...
                            Sref_Osvd, ...
                            Sref_Psvd, ...
                            Sref_Qsvd, ...
                            Sref_Gsvd];
                            % Sref_OPsvd, ...
                            % Sref_OGsvd];

Sref_svd_glb_basic_array = LD_aux.make_global_svd_package(Sref_svd_basic_array);
Sref_Lsvd_glb = Sref_svd_glb_basic_array(1);
Sref_Rsvd_glb = Sref_svd_glb_basic_array(2);
Sref_Osvd_glb = Sref_svd_glb_basic_array(3);
Sref_Psvd_glb = Sref_svd_glb_basic_array(4);
Sref_Qsvd_glb = Sref_svd_glb_basic_array(5);
Sref_Gsvd_glb = Sref_svd_glb_basic_array(6);
% Sref_OPsvd_glb = Sref_svd_glb_basic_array(7);
% Sref_OGsvd_glb = Sref_svd_glb_basic_array(8);

[Sref_OYLglbsvd,Sref_OYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_Osvd,Sref_Lsvd_glb,1);
% [Sref_OYLglbsvd,Sref_OYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_Osvd,Sref_Lsvd_glb,4);
% [Sref_OYLglbsvd,Sref_OYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_Osvd,Sref_Lsvd_glb, ...
%                                                                                                 Sref_Lsvd_glb.ncol-100);
[Sref_GYLglbsvd,Sref_GYLglbsvd_glb] = LD_observations_set.make_global_restricted_svd_package(Sref_Gsvd,Sref_Lsvd_glb,1);

% Sref_svd_Lglb_reg_array = [ Sref_OYLglbsvd ];
% Sref_svd_glb_reg_array = [ Sref_OYLglbsvd_glb ];
Sref_svd_Lglb_reg_array = [ Sref_OYLglbsvd,Sref_GYLglbsvd ];
Sref_svd_glb_reg_array = [ Sref_OYLglbsvd_glb,Sref_GYLglbsvd_glb ];

Sref_svd_glb_array = [Sref_svd_glb_basic_array, Sref_svd_glb_reg_array];

Sref_glb_svd_plot = LD_plots.plot_global_svds(Sref_svd_glb_array,LD_plots('svd',[7 7],[2 7],[5 1],1));
Sref_glb_svd_plot.show_toolbar();

[Sref_RYLsvd_in,Sref_RYLsvd_rst_in] = Sref.read_mat_svd_package(fam_name,bor,'Rmat','YL','Lmat');
Sref_svd_YL_in_array = [    Sref_RYLsvd_in];

[Sref_OYLsvd,Sref_OYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Osvd,Sref_Lsvd);
[Sref_PYLsvd,Sref_PYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Psvd,Sref_Lsvd);
[Sref_GYLsvd,Sref_GYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_Gsvd,Sref_Lsvd);
% [Sref_OGYLsvd,Sref_OGYLsvd_rst] = LD_observations_set.make_restricted_svd_package(Sref_OGsvd,Sref_Lsvd);

Sref_svd_YL_array = [       Sref_OYLsvd, ...
                            Sref_PYLsvd, ...
                            Sref_GYLsvd];
                            % Sref_OGYLsvd];

% Sref_OregRsvd = Sref.read_mat_svd_package(fam_name,bor,'OregRmat');
% Sref_svd_Oreg_array = [ Sref_OregRsvd   ];

% Sref_svd_array = [Sref_svd_basic_array, Sref_svd_Lglb_reg_array, Sref_svd_YL_in_array, Sref_svd_YL_array, Sref_svd_Oreg_array];
Sref_svd_array = [Sref_svd_basic_array, Sref_svd_Lglb_reg_array, Sref_svd_YL_in_array, Sref_svd_YL_array];
Sref_crv_svd_plot = LD_plots.plot_curve_svds(Sref_svd_array,LD_plots('svd',[7 7],[2 7],[7 1],1));
% Sref_crv_svd_plot = LD_plots.plot_curve_svds([Sref_Osvd,Sref_svd_Lglb_reg_array],LD_plots('svd',[7 7],[2 7],[6 1],1));


%%%%%%%%%%%%%%%%%%%%%%%%

Lsvd_glb = Sref_Lsvd_glb;
Asvd = Sref_Osvd;
AYLsvd = Sref_OYLsvd;
AYLrstsvd = Sref_OYLsvd_rst;
AYLglbsvd = Sref_OYLglbsvd_glb;

% Lsvd_glb = Sref_Lsvd_glb;
% Asvd = Sref_Gsvd;
% AYLsvd = Sref_GYLsvd;
% AYLrstsvd = Sref_GYLsvd_rst;
% AYLglbsvd = Sref_GYLglbsvd_glb;

[ndep,nvar] = deal(meta0.ndep,meta0.ndep+1);
[ncrv,nrow,ncol] = deal(length(Asvd.rvec),Asvd.nrow,Asvd.ncol);
ncol_glb = AYLglbsvd.ncol;
kappa_L = (ncol-ncol_glb)/nvar;
rho_L = Lsvd_glb.ncol - kappa_L;

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
% rho_PL = ncol_AYL-2;
kappa_AYL = ncol_AYL - rho_AYL
Ktns_A_YLrst = LD_aux.Ytns_Vtns_mult(Ytns_L,AYLrstsvd.Vtns(:,(rho_AYL+1):end,:)); % AYLmat kernals of each curve
Kmat_Ai_YLrst = Ktns_A_YLrst(:,:,icrv_centre); % AYLmat kernal of 'centre' curve

Kmat_Ai_T_Ktns_A = pagemtimes(Kmat_Ai_YLrst',Ktns_A_YLrst); % linear combinations of 'centre' AYLmat kernal to fit each curve's kernal
Ktns_A_K_Ai_rstr = pagemtimes(Kmat_Ai_YLrst,Kmat_Ai_T_Ktns_A); % image of each curve's kernal in terms of central kernal
mags_Ktns_A_K_Ai_rstr = reshape(sqrt(sum(Ktns_A_K_Ai_rstr.*Ktns_A_K_Ai_rstr,1)),[],ncrv); % magnitude of each image kernal

mat_colstats = @(mat_) deal(min(mat_,[],1),median(mat_,1),mean(mat_,1),max(mat_,[],1));
[min_mags,med_mags,avg_mags,max_mags] = mat_colstats(mags_Ktns_A_K_Ai_rstr);

% mag_evl = med_mags;
mag_evl = avg_mags;

[~,isort_mag_evl] = sort(mag_evl);

% icrv_neighbor = isort_mag_evl(end)

icrv_neighbor = isort_mag_evl(end-1)
mag_evl_neighbor = mag_evl(icrv_neighbor)

icrv_distant = isort_mag_evl(1)
mag_evl_distant = mag_evl(icrv_distant)


Kmat_A_YLrst = reshape(Ktns_A_YLrst,ncol,[]);
% KTK_A_YLrst = mat2cell(Kmat_A_YLrst'*Kmat_A_YLrst,kappa_AYL*ones(1,ncrv), kappa_AYL*ones(1,ncrv));
% K_A_YLrst = cell(ncrv);
% for i = 1:ncrv
%     for j = 1:ncrv
%         K_A_YLrst{i,j} = Ktns_A_YLrst(:,:,i)*KTK_A_YLrst{i,j};
%     end
% end
K_A_YLrst = mat2cell(Kmat_A_YLrst'*Kmat_A_YLrst,kappa_AYL*ones(1,ncrv), kappa_AYL*ones(1,ncrv));

mflag_lw = logical(tril(ones(ncrv),-1)); % flags of lower triangular elements
vflag_lw = reshape(mflag_lw,[],1); % column stack of lower triangular flags
vinds_full = (1:(ncrv*ncrv))'; % 1 thru ncrv^2
vinds_full_mat = reshape(vinds_full,ncrv,ncrv); % matrix of column stacked indices
vinds_lw = vinds_full(vflag_lw); % column stacked indices of strictly lower triangular elements
K_A_YLrst_lwvec = K_A_YLrst(vflag_lw); % column stacked sub cells for lower triangular elements

% avg_colmag = @(mat_) mean(sqrt(sum(mat_.*mat_,1)),2);
avg_colmag = @(mat_) sum( sqrt(sum(mat_.*mat_,1)), 2 )/( size(mat_,2)-1 );
avg_rowmag = @(mat_) avg_colmag(mat_');
% fro_nrmmag = @(mat_) norm(mat_,'fro')/sqrt(size(mat_,2));
fro_nrmmag = @(mat_) norm(mat_,'fro')/sqrt(size(mat_,2)-1);
fro2_nrmmag = @(mat_) sum(mat_(:).*mat_(:))/size(mat_,2);
nuke_nrmmag = @(mat_) mean(svd(mat_));
% nuke_nrmmag = @(mat_) sum(svd(mat_))/(size(mat_,2) - 1);
% ell2_nrmmag = @(mat_) norm(mat_);
ell2_nrmmag = @(mat_) sum(svd(mat_))/(size(mat_,2) - 1);
mins_nrmmag = @(mat_) min(svd(mat_));

avg_cmag_K_A_YLrst_lwvec = cellfun(avg_colmag,K_A_YLrst_lwvec);
avg_rmag_K_A_YLrst_lwvec = cellfun(avg_rowmag,K_A_YLrst_lwvec);
avg_fmag_K_A_YLrst_lwvec = cellfun(fro_nrmmag,K_A_YLrst_lwvec);
% avg_f2mag_K_A_YLrst_lwvec = cellfun(fro2_nrmmag,K_A_YLrst_lwvec);
% avg_nmag_K_A_YLrst_lwvec = cellfun(nuke_nrmmag,K_A_YLrst_lwvec);
% ell2_nrm_K_A_YLrst_lwvec = cellfun(ell2_nrmmag,K_A_YLrst_lwvec);

% evl_mag_K_A_YLrst_lwvec = avg_cmag_K_A_YLrst_lwvec;
evl_mag_K_A_YLrst_lwvec = avg_fmag_K_A_YLrst_lwvec;

[~,isrt_evl_mag_K_A_YLrst_lwvec] = sort(evl_mag_K_A_YLrst_lwvec); % sorted point stat of cokernal magnitudes

[idst_evl_lwvec,ingh_evl_lwvec] = deal(isrt_evl_mag_K_A_YLrst_lwvec(1),isrt_evl_mag_K_A_YLrst_lwvec(end))
[idst_evl,ingh_evl] = deal(vinds_lw(idst_evl_lwvec),vinds_lw(ingh_evl_lwvec))
[crv_dst1,crv_dst2] = find(vinds_full_mat==idst_evl)
[crv_ngh1,crv_ngh2] = find(vinds_full_mat==ingh_evl)
distant_evl_mag_K_A_YLrst_lwvec = evl_mag_K_A_YLrst_lwvec(idst_evl_lwvec)
neighbor_evl_mag_K_A_YLrst_lwvec = evl_mag_K_A_YLrst_lwvec(ingh_evl_lwvec)


K_A_YLrst_evl_mag_lw = zeros(ncrv);
K_A_YLrst_evl_mag_lw(vflag_lw) = evl_mag_K_A_YLrst_lwvec;
K_A_YLrst_evl_mag = K_A_YLrst_evl_mag_lw + eye(ncrv) + K_A_YLrst_evl_mag_lw';

[min_evlmags,med_evlmags,avg_evlmags,max_evlmags] = mat_colstats(K_A_YLrst_evl_mag);

evl_evlmags = avg_evlmags;
% evl_evlmags = med_evlmags;
% evl_evlmags = max_evlmags;

[~,isrt_evl_evlmags] = sort(evl_evlmags);
% [icrv_loneliest,icrv_frndliest] = deal(isrt_evl_evlmags(1),isrt_evl_evlmags(end))
[icrv_loneliest,icrv_frndliest] = deal(isrt_evl_evlmags(2),isrt_evl_evlmags(end))

[~,isort_mag_evl2] = sort(K_A_YLrst_evl_mag(:,icrv_frndliest));
icrv_neighbor2 = isort_mag_evl2(end-1)
icrv_distant2 = isort_mag_evl2(1)

%%%%%

spc.lw = 2.0;
% spc.color = [0 0 0];
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_centre);
% spc.color = [0 0 1];
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_neighbor);
% spc.color = [1 0 0];
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_distant);

spc.color = LD_plots.orange1;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,crv_dst1);
spc.color = LD_plots.orange2;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,crv_dst2);

spc.color = LD_plots.blue1;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,crv_ngh1);
spc.color = LD_plots.blue2;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,crv_ngh2);

spc.color = LD_plots.purple1;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_loneliest);
spc.color = LD_plots.purple5;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_frndliest);

% spc.color = LD_plots.blue4;
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_neighbor2);
% spc.color = LD_plots.red3;
% solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_distant2);

solspc_plot_init_lims = solspc_ref_plot.get_axis_lims;

aux_fig1 = LD_plots('aux1',[5 5],[1 5],[5 1],1);

aux_fig1 = aux_fig1.init_tiles_safe(1,4);
hold(aux_fig1.axs, 'on');
box(aux_fig1.axs,'on');
axs_aux1 = aux_fig1.axs;

[~,isrt_pltmg] = sort(avg_fmag_K_A_YLrst_lwvec);

plot(axs_aux1(1),1:size(avg_cmag_K_A_YLrst_lwvec,1),avg_cmag_K_A_YLrst_lwvec(isrt_pltmg), ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 0 1]);
plot(axs_aux1(1),1:size(avg_rmag_K_A_YLrst_lwvec,1),avg_rmag_K_A_YLrst_lwvec(isrt_pltmg), ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[1 0 0]);
plot(axs_aux1(1),1:size(avg_fmag_K_A_YLrst_lwvec,1),avg_fmag_K_A_YLrst_lwvec(isrt_pltmg), ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 1 0]);
% plot(axs_aux1(1),1:size(avg_f2mag_K_A_YLrst_lwvec,1),avg_f2mag_K_A_YLrst_lwvec(isrt_pltmg), ...
% 'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',LD_plots.green4);
% plot(axs_aux1(1),1:size(avg_nmag_K_A_YLrst_lwvec,1),avg_nmag_K_A_YLrst_lwvec(isrt_pltmg), ...
% 'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',LD_plots.blue1);
% plot(axs_aux1(1),1:size(ell2_nrm_K_A_YLrst_lwvec,1),ell2_nrm_K_A_YLrst_lwvec(isrt_pltmg), ...
% 'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',LD_plots.orange1);

plot(axs_aux1(2),avg_fmag_K_A_YLrst_lwvec,avg_cmag_K_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 0 1]);
plot(axs_aux1(2),avg_fmag_K_A_YLrst_lwvec,avg_rmag_K_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[1 0 0]);

plot(axs_aux1(3),avg_rmag_K_A_YLrst_lwvec,avg_cmag_K_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 0 1]);
plot(axs_aux1(3),avg_rmag_K_A_YLrst_lwvec,avg_fmag_K_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 1 0]);

plot(axs_aux1(4),avg_cmag_K_A_YLrst_lwvec,avg_rmag_K_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[1 0 0]);
plot(axs_aux1(4),avg_cmag_K_A_YLrst_lwvec,avg_fmag_K_A_YLrst_lwvec, ...
'Marker','o','MarkerSize',5,'LineStyle','none','LineWidth',1,'Color',[0 1 0]);

aux_fig1.show_toolbar();

return

%%%%%%%%%%%%%%%%%%%%%%%%

% Kglb = AYLglbsvd.Vmat_glb(:,end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-ceil(ncol_glb/2)):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-50):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-40):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-30):end);
Kglb = AYLglbsvd.Vmat_glb(:,(end-24):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-22):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-20):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-16):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-14):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-12):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-10):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-8):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-7):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-6):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-4):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-2):end);
% Kglb = AYLglbsvd.Vmat_glb(:,(end-1):end);

kappa_Kglb = size(Kglb,2)

YLKglb = LD_aux.Ytns_Vtns_mult(Lsvd_glb.Vmat_glb(:,1:rho_L),Kglb);
absAtnsYLKglb = abs(pagemtimes(Atns,YLKglb));
absinn_AYLKg = reshape(sum(absAtnsYLKglb,1),size(absAtnsYLKglb,2),[]);

absinn_AYLKg_evl = sum(absinn_AYLKg,1);
% absinn_AYLKg_evl = median(absinn_AYLKg,1);

[~,isort_absinn_AYLKg_evl] = sort(absinn_AYLKg_evl);

icrv_check = isort_absinn_AYLKg_evl(1)

evl_icrv_check = absinn_AYLKg_evl(icrv_check)

pts_cell = Sref.pts_cell;
pts_mat_i = pts_cell{icrv_check};

vcurve_check = fspace0.ds_de_ratio_ptsmat_stabilized_fast(pts_mat_i,YLKglb);
scurve_check = [pts_mat_i(1:nvar,:) ; vcurve_check(2:(end-ndep),:)];

spc.color = [0 0 0];
spc.lw = 2.0;
solspc_ref_plot = LD_plots.plot_solspc(Sref,solspc_ref_plot,spc,icrv_check);

solspc_plot_init_lims = solspc_ref_plot.get_axis_lims;

spc.color = [0 0 1];
solspc_plot_ = LD_plots.plot_pts(scurve_check,meta0,solspc_ref_plot,spc);

solspc_ref_plot.set_axis_lims(solspc_plot_init_lims);
