classdef LD_denoise_plots < LD_plots
    properties

    end
    methods
        function obj = LD_denoise_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_)
            obj@LD_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_);
        end
    end
    methods (Static)
        function plt_out = plot_denoised_trajectories(plt_,spc_,Sref_,Snse_,icrv_)
            if (nargin==4)
                icrv_plot = 1:Sref_.ncrv;
            else
                icrv_plot = icrv_;
            end
            plt_out = plt_;
            spc = spc_;

            meta0 = Sref_.meta_data();

            % icrv_plot = 1:Sref_.ncrv;

            ndns_short_range = 1:5;

            % ndns_range = 1:3;
            % ndns_range = [ 1:5 ];
            % ndns_range = [ 1:5 , 6:2:8 ];
            % ndns_range = [ 1:5 , 6:2:10 ];
            % ndns_range = [ 1:5 , 6:2:40 ];
            % ndns_range = [ 1:5 , 6:2:50 ];
            % ndns_range = [ 1:5 , 6:2:100 ];
            % ndns_range = [ 1:5 , 10:10:100 ];
            % ndns_range = [ 1:5 , 10:10:200 ];
            % ndns_range = [ 1:5 , 10:10:300 ];
            % ndns_range = [ 1:5 , 10:10:400 ];
            % ndns_range = [ 1:5 , 10:10:90, 94 ];
            % ndns_range = [ 1:5 , 10:10:90, 99 ];

            % ndns_range = [ 1:5 , 20:20:500 ];
            % ndns_range = [ 1:5 , 50:50:600 ];
            % ndns_range = [ 1:5 , 100:100:900 , 999];
            % ndns_range = [ 1:5 , 100:100:200 , 210];
            % ndns_range = [ 1:5 , 100:100:200 , 210];

            ndns_long_range = (ndns_short_range(end)+1):1:10;

            ndns_range = [ndns_short_range ndns_long_range];

            len_ndns_range = length(ndns_range);

            pts_cell = cell([Sref_.ncrv,len_ndns_range]);

            color_mat = flip(cool(len_ndns_range),1);

            color_mat_short = flip(spring( length(ndns_short_range) ),1);
            color_mat_long = flip(cool( length(ndns_long_range) ),1);
            color_mat = [color_mat_short ; color_mat_long];

            [spc.lspec,spc.lw] = deal('-',0.5);
            [spc.mspec,spc.ms] = deal('o',3);
            for i = 1:len_ndns_range
                idns = ndns_range(i);
                pts_cell(:,i) = Snse_.read_Sobs_cell(['.jsol_Rk_' num2str(idns)]);

            % spc.color = [color_mat(i,:) 0.5];
            % plt_out = LD_plots.plot_pts(Snse_.read_Sobs_cell(['.jsol_R1_' num2str(idns)]), ...
            %                                 meta0,plt_out,spc);

            end

            % [spc.lspec,spc.lw] = deal('-',1);
            % [spc.mspec,spc.ms] = deal('.',1);
            % spc.color = [0 0 0 0.25];
            % plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc);

                % [spc.lspec,spc.lw] = deal('none',0.5);
                % [spc.mspec,spc.ms] = deal('s',3);
                % spc.color = [1 0 0 0.25];
                % plt_out = LD_plots.plot_solspc(Snse_,plt_out,spc);

            [spc.lspec,spc.lw] = deal('-',0.5);
            [spc.mspec,spc.ms] = deal('.',5);
            spc.color = [0 0 0 1];
            plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc,icrv_plot);

            [spc.lspec,spc.lw] = deal('-',0.5);
            [spc.mspec,spc.ms] = deal('s',3);
            spc.color = [1 0 0 1];
            plt_out = LD_plots.plot_solspc(Snse_,plt_out,spc,icrv_plot);


        % axlim = plt_out.get_axis_lims();
                % [spc.lspec,spc.lw] = deal('-',0.5);
                % [spc.mspec,spc.ms] = deal('none',5);
                % spc.color = [0 0 0 0.5];
                % plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc);
        % plt_out.set_axis_lims(axlim);

            icrv_plot = 1:6;
            % icrv_plot = 1:10;
            % icrv_plot = 1:50;

            [spc.lspec,spc.lw] = deal('-',1);
            [spc.mspec,spc.ms] = deal('.',1);
            spc.color = [0 0 0 0.25];
            plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc,icrv_plot);

            [spc.lspec,spc.lw] = deal('none',0.5);
            [spc.mspec,spc.ms] = deal('s',3);
            spc.color = [1 0 0 0.25];
            plt_out = LD_plots.plot_solspc(Snse_,plt_out,spc,icrv_plot);

            [spc.lspec,spc.lw] = deal('-',0.5);
            [spc.mspec,spc.ms] = deal('o',3);

            for i = 1:(len_ndns_range-1)
                spc.color = [color_mat(i,:) 1];
                plt_out = LD_plots.plot_pts(pts_cell(icrv_plot,i), ...
                                        meta0,plt_out,spc);

            end

            [spc.lspec,spc.lw] = deal('-',0.5);
            [spc.mspec,spc.ms] = deal('*',3);
            spc.color = [color_mat(len_ndns_range,:) 1];
            plt_out = LD_plots.plot_pts(pts_cell(icrv_plot,len_ndns_range), ...
                                    meta0,plt_out,spc);


            % [spc.mspec,spc.ms] = deal('s',4);
            % [spc.lspec,spc.lw] = deal('none',1);
            % % spc.color = [1 0 0 1];
            % spc.color = LD_plots.orange1;
            % plt_out = LD_plots.plot_solspc(Snse_,plt_out,spc,icrv_plot);
        end
        function plt = plot_rowspace_image(plt_,spc_,Sref_,Snse_,svd_,rmg_,svd_n_,rmg_n_,rmg_st_n_)
            plt = plt_;
            [tdim1,tdim2] = deal(3,2);
            plt = plt.init_tiles_safe(tdim1,tdim2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = plt.axs_mat();

            pts_mat = Sref_.pts_mat;
            pts_mat_n = Snse_.pts_mat;

            x_vec = pts_mat(1,:);

            dnp1xu_mat = Sref_.dnp1xu_mat();
            dnp1xu_mat_n = Snse_.dnp1xu_mat();

            vnu_mat = rmg_.vnumat;
            vnu_mat_n = rmg_n_.vnumat;
            vnu_mat_st_n = rmg_st_n_.vnumat;

            vnu_s_mat = rmg_.vnusmat;
            vnu_s_mat_n = rmg_n_.vnusmat;
            vnu_s_mat_st_n = rmg_st_n_.vnusmat;

            ichk = 1:9;
            snp1_chk = [pts_mat(:,ichk) ; dnp1xu_mat(:,ichk)];
            vnu_tru = [pts_mat((Sref_.ndim-Sref_.ndep+1):end,ichk) ; dnp1xu_mat(:,ichk)]

            vnu_chk = vnu_mat(:,ichk)
            vnu_s_chk = vnu_s_mat(:,ichk)

            vnu_st_n_chk = vnu_mat_st_n(:,ichk)
            vnu_s_st_n_chk = vnu_s_mat_st_n(:,ichk)

            snp1_n_chk = [pts_mat_n(:,ichk) ; dnp1xu_mat_n(:,ichk)];
            vnu_nse = [pts_mat_n((Sref_.ndim-Sref_.ndep+1):end,ichk) ; dnp1xu_mat_n(:,ichk)]

            vnu_tru_full = [pts_mat((Sref_.ndim-Sref_.ndep+1):end,:) ; dnp1xu_mat];
            vnu_nse_full = [pts_mat_n((Sref_.ndim-Sref_.ndep+1):end,:) ; dnp1xu_mat_n];
            vnu_ref_full = vnu_mat;
            vnu_ref_s_full = vnu_s_mat;
            vnu_nsn_full = vnu_mat_n;
            vnu_nsn_s_full = vnu_s_mat_n;
            vnu_nst_full = vnu_mat_st_n;
            vnu_nst_s_full = vnu_s_mat_st_n;

            err = @(y_,k_) y_(k_,:) - vnu_tru_full(k_,:);
            err_nse = @(y_,k_) y_(k_,:) - vnu_nse_full(k_,:);
            abserr = @(y_,k_) abs(err(y_,k_));

            exp_plot = {
                        vnu_nst_full, LD_plots.blue5 ; ...
                        vnu_nsn_full, LD_plots.purple5 ; ...
                        vnu_nse_full, LD_plots.red5 ; ...
                        vnu_ref_full, LD_plots.green4 ; ...
                       };

            spc.mspec = '.';
            spc.lspec = 'none';
            spc.ms = 6;
            spc.lw = 0.5;
            spc.alpha = 0.25;

            axi = axs(1);

            nexp = size(exp_plot,1);

            k = 1;
            for i = 1:nexp
                plot(axs(1), ...
                    x_vec, ...
                    abserr(exp_plot{i,1},k), ...
                    'Color',exp_plot{i,2}, ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
                % scatter(axs(1), ...
                %     x_vec, ...
                %     abserr(exp_plot{i,1},k), ...
                %     'Marker', spc.mspec, ...
                %     'MarkerFaceColor',exp_plot{i,2}, ...
                %     'MarkerFaceAlpha',spc.alpha, ...
                %     'MarkerEdgeAlpha',spc.alpha, ...
                %     'SizeData',spc.ms,'LineWidth',spc.lw)
            end
            % for i = 1:nexp
            % for i = 2:nexp
            % for i = 3:nexp
            for i = 1:(nexp-1)
                histogram(axs(2), ...
                    err(exp_plot{i,1},k), ...
                    'Normalization', 'probability', ...
                    'FaceColor',exp_plot{i,2}, 'FaceAlpha',0.4);
            end

            for i = 1:nexp
                plot(axs(3), ...
                    x_vec, ...
                    err(exp_plot{i,1},k), ...
                    'Color',exp_plot{i,2}, ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end
            for i = 1:nexp
                plot(axs(4), ...
                    x_vec, ...
                    err_nse(exp_plot{i,1},k), ...
                    'Color',exp_plot{i,2}, ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end
            for i = 1:nexp
                plot(axs(5), ...
                    err_nse(exp_plot{i,1},k), ...
                    err(exp_plot{i,1},k), ...
                    'Color',exp_plot{i,2}, ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end


            set(axs, 'TickLabelInterpreter','Latex','FontSize',12);
            set(axs(1),'YScale','log', 'XScale','linear');
            % set(axs(5),'YScale','log', 'XScale','log');

        end
        function plt = plot_Amat_data(plt_,spc_,Aref_,Asvd_ref_,Anse_,Asvd_nse_)
            plt = plt_;
            [tdim1,tdim2] = deal(1,3);
            plt = plt.init_tiles_safe(tdim1,tdim2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = (reshape(axs,tdim2,tdim1))';

            Smat_ref = cell2mat(Asvd_ref_(5,:))';
            Smat_nse = cell2mat(Asvd_nse_(5,:))';
            [nspc,vlen] = size(Smat_ref);

            spc = spc_;

            [AV_ref,AV_nse,ArVn,AnVr] = deal(nan(nspc,vlen));
            for i = 1:nspc
                AV_ref(i,:) = sqrt(sum(((Aref_(i).mat)*Asvd_ref_{4,i}).^2,1));
                AV_nse(i,:) = sqrt(sum(((Anse_(i).mat)*Asvd_nse_{4,i}).^2,1));
                ArVn(i,:) = sqrt(sum(((Aref_(i).mat)*Asvd_nse_{4,i}).^2,1));
                AnVr(i,:) = sqrt(sum(((Anse_(i).mat)*Asvd_ref_{4,i}).^2,1));
            end

            plt_cell = {Smat_ref,Smat_nse ; AV_ref,AV_nse ; ArVn,AnVr };
            color_mat = [LD_plots.blue5 ; LD_plots.red5];
            for iax = 1:size(plt_cell,1)
                for i = 1:size(plt_cell,2)
                    plot(axs(iax), ...
                        1:vlen,median(plt_cell{iax,i},1), ...
                        '-','Color',color_mat(i,:),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                    plot(axs(iax), ...
                        1:vlen,min(plt_cell{iax,i},[],1), ...
                        ':','Color',color_mat(i,:),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                    plot(axs(iax), ...
                        1:vlen,max(plt_cell{iax,i},[],1), ...
                        '--','Color',color_mat(i,:),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                end
            end

            % plot(axs(3), ...
            %     Smat_ref(:),AV_ref(:), ...
            %     'o','Color',color_mat(1,:),'MarkerSize',spc.ms,'LineWidth',spc.lw)
            % plot(axs(3), ...
            %     Smat_nse(:),AV_nse(:), ...
            %     's','Color',color_mat(2,:),'MarkerSize',spc.ms,'LineWidth',spc.lw)


            set(axs(:), 'TickLabelInterpreter','Latex','FontSize',12);
            % set(axs(1:3),'YScale','log', 'XScale','linear');
            set(axs(1:3),'YScale','linear', 'XScale','linear');
            % set(axs(3),'YScale','log', 'XScale','log');

            % LD_plots.set_containing_axis_lims(axs_mat(1,1))
            % LD_plots.set_containing_axis_lims(axs_mat(2,1))

        end
        function plt = plot_curve_estimates(Splt_,plt_,spc_,Sref_,Snse_,icrv_)
            plt = plt_;
            [tdim1,tdim2] = deal(3,Sref_.eor+1);
            plt = plt.init_tiles_safe(tdim1,tdim2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = (reshape(axs,tdim2,tdim1))';

            meta0 = LD_observations_set.make_meta_data(Sref_.eor,Sref_.ndep);
            spc = spc_;
            [eor,ndep] = deal(Sref_.eor,Sref_.ndep);
            i_crv = icrv_;
            kor_u_hat = eor;

            sigma_un = reshape(Snse_.sigma_s(2:end),ndep,[])';
            sigma_unp1 = [  sigma_un; ...
                            reshape(Snse_.sigma_dnp1xu,1,ndep) ];

            lmat_kp1 = ones(ndep,1);
            lmat_crv = ones(eor+1,ndep);

            % sigma_scl = 1e-8;
            % sigma_scl = 1e-7;
            % sigma_scl = 1e-6;
            % sigma_scl = 1e-5;
            % sigma_scl = 1e-4;
            % sigma_scl = 1e-3;
            % sigma_scl = 1e-2;
            % sigma_scl = 1e-1;
            sigma_scl = 1e0;

            % lam_scl = 1e-8;
            % lam_scl = 1e-6;
            lam_scl = 1e0;
            % lam_scl = 1e1;
            % lam_scl = 1e2;
            % lam_scl = 1e3;
            % lam_scl = 1e4;
            % lam_scl = 1e6;
            % lam_scl = 1e8;
            % lam_scl = 1e10;

            s_un_chk = sigma_scl*sigma_un;
            s_unp1_chk = sigma_scl*sigma_unp1;
            lmat_kp1_chk = lam_scl*lmat_kp1;
            lmat_crv_chk = lam_scl*lmat_crv;

            [crvs_ref_,crvs_nse_] = deal(Sref_.crvs,Snse_.crvs);

            crv_i = crvs_ref_(i_crv);
            jet_i = crv_i.jets;
            jet_in = crv_i.make_jets(kor_u_hat);

            crv_i_n = crvs_nse_(i_crv);
            jet_i_n = crv_i_n.jets;
            jet_in_n = crv_i_n.make_jets(kor_u_hat);

            rsjet_i = LD_jets.regularize_jet(jet_i,s_unp1_chk,lmat_kp1_chk);
            rcjet_i = LD_jets.minavg_curvature_jet(jet_i,s_unp1_chk,lmat_crv_chk);
            rsjet_i_n = LD_jets.regularize_jet(jet_i_n,s_unp1_chk,lmat_kp1_chk);
            rcjet_i_n = LD_jets.minavg_curvature_jet(jet_i_n,s_unp1_chk,lmat_crv_chk);

            rsjet_in = LD_jets.regularize_jet(jet_in,s_un_chk,lmat_kp1_chk);
            rcjet_in = LD_jets.minavg_curvature_jet(jet_in,s_un_chk,lmat_crv_chk);
            rsjet_in_n = LD_jets.regularize_jet(jet_in_n,s_un_chk,lmat_kp1_chk);
            rcjet_in_n = LD_jets.minavg_curvature_jet(jet_in_n,s_un_chk,lmat_crv_chk);

            % x_check = 0:0.01:25;
            % x_check = jet_i.pts_mat(1,:);
            % x_check = jet_i.xh_vec;
            x_check = sort([jet_i.xh_vec,jet_i.pts_mat(1,:)]);
            inds_01 = 1:(2):length(x_check);
            inds_h = 2:(2):(length(x_check)-1);

            jexps_full = { jet_i, [] , LD_plots.green5 ; ... % (1) unnoised Hermite jets, kor = eor+1
                    jet_in, [], [0 0 0] ; ... % (2) unnoised Hermite jets, kor = eor
                jet_i_n, [], LD_plots.red5 ; ... % (3) noisy Hermite jets, kor = eor+1
                    jet_in_n, [], LD_plots.orange1 ; ... % (4) noisy Hermite jets, kor = eor
                rsjet_i, [], LD_plots.blue5; ... % (5) unnoised kor+1 jets, kor = eor+1
                    rsjet_in, [], LD_plots.blue1; ... % (6) unnoised kor+1 jets, kor = eor
                rsjet_i_n, [], LD_plots.orange3; ... % (7) noisy kor+1 reg. jets, kor = eor+1
                    rsjet_in_n, [], LD_plots.orange5; ... % (8) noisy kor+1 reg. jets, kor = eor
                rcjet_i, [], LD_plots.green4 ; ... % (9) unnoised crv. reg. jets, kor = eor+1
                    rcjet_in, [], LD_plots.green1 ; ... % (10) unnoised crv. reg. jets, kor = eor
                rcjet_i_n, [], LD_plots.purple1 ; ... % (11) noisy crv. reg. jets, kor = eor+1
                    rcjet_in_n, [], LD_plots.purple5 ; ... % (12) noisy crv. reg. jets, kor = eor
                        };
            njets_full = size(jexps_full,1);

            for i = 1:njets_full
                [u_check_i,dxu_check_i] = jexps_full{i,1}.u_hat(x_check,eor);
                jexps_full{i,2} = [u_check_i ; dxu_check_i];
            end
            Uref = jexps_full{1,2};

            jexps_plot = jexps_full([   1, ... % jexp_i
                                        3, ... % jexp_i_n
                                        7, ... % rsjexp_i_n
                                        11 ... % rcjexp_i_n
                                    ], :);
            % jexps_plot = jexps_full([   2, ... % jexp_in
            %                             8, ... % rsjexp_in_n
            %                             12 ... % rcjexp_in_n
            %                         ], :);

            njexp_plot = size(jexps_plot,1);

            spc.mspec = '.';
            spc.ms = 8; % spc.ms = 8;
            for i = 1:njexp_plot
                spc.color = jexps_plot{i,3};
                LD_plots.plot_pts([x_check; jexps_plot{i,2}],meta0,Splt_,spc);
            end

            zscore = @(u_,y_,k_) (y_-u_)/sigma_unp1(k_);
            zscabs = @(u_,y_,k_) abs(zscore(u_,y_,k_));
            zsqurd = @(u_,y_,k_) zscore(u_,y_,k_).^2;
            rawres = @(u_,y_) y_-u_;
            absres = @(u_,y_) abs(u_-y_);
            abserr = @(u_,y_) abs((u_-y_)./u_);
            mdplt = @(x_) ones(1,2)*median(x_);
            medabserr = @(u_,y_) abs((u_-y_)./u_);
            row_stats = @(row_) deal(min(row_,[],2),median(row_,2),mean(row_,2),max(row_,[],2));
            udkstr = @(k_) ['u^{[' num2str(k_-1) ']}'];
            err_str = @(k_) ['\psi^{' udkstr(k_) '}-' udkstr(k_)];
            zscore_str = @(k_) [ '(' err_str(k_) ')/ \sigma_{' udkstr(k_) '}'];

            icll = {    inds_01,'v','--'; ...
                        inds_h,'^',':'; ...
                    };
            ninds = size(icll,1);
            % evl = { @(u_,y_,k_) absres(u_,y_) ; ...
            evl = { @(u_,y_,k_) zscabs(u_,y_,k_), ...
                     @(k_) ['$$ |' zscore_str(k_) '| $$']; ...
                    @(u_,y_,k_) zscore(u_,y_,k_), ...
                     @(k_) ['$$' zscore_str(k_) '$$']; ...
                    };

            % jexps_eval = jexps_plot;
            jexps_eval = jexps_full([   3, ... % jet_i_n
                                        7, ... % rsjexp_i_n
                                        11, ... % rcjexp_i_n
                                    ], :);
            % jexps_eval = jexps_full([   4, ... % jet_in_n
            %                             8, ... % rsjexp_in_n
            %                             12, ... % rcjexp_in_n
            %                         ], :);
            nexp_eval = size(jexps_eval,1);

            spc.mspec = spc_.mspec;
            spc.ms = 2; % spc.ms = spc_.ms;
            for ievl = 1:size(evl,1)
                for ik = 1:(eor+1)
                    axi = axs_mat(ievl,ik);
                    for iexp = 1:nexp_eval
                        for jinds = 1:ninds
                            if (~( (iexp==1)&&(jinds==2)&&(ik>1) ))
                                plot(axi, ...
                                    x_check(icll{jinds,1}), ...
                                    evl{ievl,1}(Uref(ik,icll{jinds,1}), ...
                                        jexps_eval{iexp,2}(ik,icll{jinds,1}),ik), ...
                                    'Marker', icll{jinds,2}, ...
                                    'LineStyle', 'none', ...
                                    'Color',jexps_eval{iexp,3}, ...
                                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
                                plot(axi, ...
                                    [x_check(1),x_check(end)], ...
                                    mdplt(evl{ievl,1}(Uref(ik,icll{jinds,1}), ...
                                        jexps_eval{iexp,2}(ik,icll{jinds,1}),ik) ), ...
                                    'Marker', 'none', ...
                                    'LineStyle', icll{jinds,3}, ...
                                    'Color',jexps_eval{iexp,3}, ...
                                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
                            end
                        end
                    end
                end
            end
                hstll = cell([eor+1,nexp_eval,ninds]);
                for ik = 1:(eor+1)
                    axi = axs_mat(3,ik);
                    for iexp = 1:nexp_eval
                        for jinds = 1:ninds
                            if (~( (iexp==1)&&(jinds==2)&&(ik>1) ))
                                hstll{ik,iexp,jinds} = histogram(axi, ...
                                    zscore(Uref(ik,icll{jinds}), ...
                                    jexps_eval{iexp,2}(ik,icll{jinds,1}),ik), ...
                                    'Normalization', 'probability', ...
                                    'FaceColor',jexps_eval{iexp,3}, 'FaceAlpha',0.5);
                            end
                            % pause
                        end
                    end
                end
                % set([hstll{2:3,1,2}],'FaceAlpha',0)
            % pause
            axs_set = axs_mat(:);
            set(axs_set,  ...
                'TickLabelInterpreter','Latex', ...
                'FontSize',12);
            set(axs_mat(1,:),'YScale','log', 'XScale','linear');
            set(axs_mat(2,:),'YScale','linear', 'XScale','linear');

            xlabel(axs_mat(2,:),'$$ x $$','Interpreter','Latex','FontSize',16)
            for k = 1:(eor+1)
                for ievl = 1:size(evl,1)
                    ylabel(axs_mat(ievl, k),evl{ievl,2}(k), ...
                        'Interpreter','Latex','FontSize',16)
                end
                xlabel(axs_mat(end, k),['$$' zscore_str(k) '$$'], ...
                        'Interpreter','Latex','FontSize',16)
            end

            % LD_plots.set_containing_axis_lims(axs_mat(1,1))
            % LD_plots.set_containing_axis_lims(axs_mat(2,1))

        end
        function plt_out = plot_observed_trajectories(plt_,spc_,Sref_,Snse_,icrv_)
            if (nargin==4)
                icrv_plot = 1;
            else
                icrv_plot = icrv_;
            end
            plt_out = plt_;
            spc = spc_;

            % spc.lw = 1;
            [spc.lspec,spc.lw] = deal('-',1);
            [spc.mspec,spc.ms] = deal('.',spc.lw);
            % plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc);

            [spc.lspec,spc.lw] = deal('none',1);
            [spc.mspec,spc.ms] = deal('o',4);
            % spc.color = [0 0 1 1];
            spc.color = [0 0 0 1];
            plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc,icrv_plot);

            [spc.mspec,spc.ms] = deal('s',4);
            [spc.lspec,spc.lw] = deal('none',1);
            % spc.color = [1 0 0 1];
            spc.color = LD_plots.orange1;
            plt_out = LD_plots.plot_solspc(Snse_,plt_out,spc,icrv_plot);
        end
    end
end
