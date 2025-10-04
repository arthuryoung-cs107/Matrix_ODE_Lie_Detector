classdef LD_denoise_plots < LD_plots
    properties

    end
    methods
        function obj = LD_denoise_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_)
            obj@LD_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_);
        end
    end
    methods (Static)
        function [plts_out,dnse_data] = plot_denoising_summary(plts_,spc_,Ssets_,ode_,icrv_)
            if (nargin==3)
                icrv_plot = 1:(Ssets_(1).ncrv);
            else
                icrv_plot = icrv_;
            end
            res_vec = [];
            plts_out = plts_;

            % Rsvd_g_names = { 'Rsvd_g' ; 'Rsvd_h_g' };
            Rsvd_g_names = { 'Rsvd_g' };
            jet_sol_names = { '.jsol_h'; ...
            '.jsol_h_Rk'; ...
            '.jsol_0_Rk'; ...
            '.jsol_1_Rk'; ...
            };

            Sref = Ssets_(1);
            meta0 = Sref.meta_data();
            [eor,ndep,ndim] = deal(meta0.eor,meta0.ndep,meta0.ndim);
            ncrv = Sref.ncrv;

            [Rsvd_g_0_tru,pSj_0_tru] = Sref.read_jet_sol_h_data(Rsvd_g_names, jet_sol_names);
            [Rsvd_g_s_tru,pSj_s_tru] = Sref.read_jet_sol_h_data(Rsvd_g_names, jet_sol_names, '_s');
            [Rsvd_g_f_tru,pSj_f_tru] = Sref.read_jet_sol_h_data(Rsvd_g_names, jet_sol_names, '_f');
            SVT_g_0_tru = (Rsvd_g_0_tru.s)' .* (Rsvd_g_0_tru.V)';

            Snse = Ssets_(2);
            % Snse = Sref;

            pts_ref_cell = Sref.pts_cell();
            pts_nse_cell = Snse.pts_cell();

            % startup_data = Snse.read_startup_data('_s');
            dnse_summary = Snse.read_denoise_summary('.denoise_summary')
            Sobs_f = Snse.read_Sobs_cell(['.jsol_Rk_' num2str(dnse_summary.nsmooth)]);

            iwrite_full = dnse_summary.iwrite;
            len_iwrite_full = length(iwrite_full);
            nwrite_short_try = 5;
            if (len_iwrite_full>=nwrite_short_try)
                iwrite_short = iwrite_full(1:nwrite_short_try);
            else
                % iwrite_short = [];
                iwrite_short = iwrite_full;
            end
            nwrite_short = length(iwrite_short);
            residuals_iwrite = dnse_summary.residuals(iwrite_full);
            ranks_iwrite = dnse_summary.ranks(iwrite_full);

            ii_iwrite_plot = logical([ zeros(len_iwrite_full-1, 1) ; 1]); % always plot last iteration
            ii_iwrite_plot(1:nwrite_short) = true; % plot first few iterations
            ii_iwrite_plot( logical([true ; diff(ranks_iwrite)<0 ]) ) = true; % plot any iterations with big errors
            ii_iwrite_plot(residuals_iwrite > 1e-1) = true; % plot any iterations with big errors
            ii_iwrite_plot( [true ; diff(residuals_iwrite) > 0.0] ) = true; % plot any increasing iterations

            inds_iwrite_full = 1:len_iwrite_full;
            iplot_iwrite = inds_iwrite_full(ii_iwrite_plot);
            iwrite_plot = iwrite_full(iplot_iwrite);

            residuals_iplot_iwrite = residuals_iwrite(iplot_iwrite);
            ranks_iplot_iwrite = ranks_iwrite(iplot_iwrite);

            len_iwrite_plot = length(iwrite_plot);
            iwrite_long = iwrite_plot((nwrite_short+1):end);
            nwrite_long = length(iwrite_long);

            [Rsvd_g_0,pSj_0] = Snse.read_jet_sol_h_data(Rsvd_g_names, jet_sol_names);
            [Rsvd_g_s,pSj_s] = Snse.read_jet_sol_h_data(Rsvd_g_names, jet_sol_names, '_s');
            [Rsvd_g_f,pSj_f] = Snse.read_jet_sol_h_data(Rsvd_g_names, jet_sol_names, '_f');
            [dxuk_Rk,dxuk_Rk_cell] = Snse.read_dxuk_data('.dxuk_Rk');
            pSj_Rk_cell = LD_observations_set.substitute_dkxu_pts_cells(pts_nse_cell,dxuk_Rk,dxuk_Rk_cell);

            ndof = length(Rsvd_g_0.s);

            Sdnse_cell = cell([Snse.ncrv,len_iwrite_plot]);
            Smat = nan(length(Rsvd_g_0.s), len_iwrite_plot);
            Vtns = nan([ size(Rsvd_g_0.V), len_iwrite_plot]);
            for i = 1:len_iwrite_plot
                Sdnse_cell(:,i) = Snse.read_Sobs_cell(['.jsol_Rk_' num2str(iwrite_plot(i))]);
                svd_i = Snse.read_LD_svd([Rsvd_g_names{1} '_' num2str(iwrite_plot(i))]);
                [Smat(:,i),Vtns(:,:,i)] = deal(svd_i.s,svd_i.V);
            end

            err_evl_tru = @(ptsS_,j_) ptsS_{j_}(2:end,:)-pts_ref_cell{j_}(2:end,:);
            err_evl_nse = @(ptsS_,j_) ptsS_{j_}(2:end,:)-pts_nse_cell{j_}(2:end,:);
            err_evl_sys = @(ptsS_,j_) ptsS_{j_}((end-ndep+1):end,:)-ode_.dnxu(ptsS_{j_}(1:(ndim-ndep),:));

            [errSnse_tru,errSnse_sys] = deal(cell([ncrv,1]));
            res_nse_ref = zeros( Snse.ndep*(Snse.eor+1),1 ) ;
            res_nse_sys = zeros( Snse.ndep, 1 ) ;
            for j = 1:ncrv
                errSnse_tru{j} = err_evl_tru( pts_nse_cell, j );
                errSnse_sys{j} = err_evl_sys( pts_nse_cell, j );

                res_nse_ref = res_nse_ref + ...
                sum( ( ...
                        pts_nse_cell{j}(2:end,:) ...
                        - pts_ref_cell{j}(2:end,:) ...
                    ).^2 ,2);
                res_nse_sys = res_nse_sys + ...
                sum( ( ...
                        pts_nse_cell{j}((end-ndep+1):end,:) ...
                        - ode_.dnxu(pts_nse_cell{j}(1:(ndim-ndep),:)) ...
                    ).^2 ,2);
            end
            err_nse_ref = sqrt(res_nse_ref);
            err_nse_sys = sqrt(res_nse_sys);

            [errSdnse_tru,errSdnse_nse,errSdnse_sys] = deal(cell([ncrv,len_iwrite_full]));
            res_tru = zeros( Snse.ndep*(Snse.eor+1),len_iwrite_full ) ;
            res_sys = zeros( Snse.ndep, len_iwrite_full ) ;
            for i = 1:len_iwrite_full
                Sdnse_cell_i = Snse.read_Sobs_cell(['.jsol_Rk_' num2str(iwrite_full(i))]);

                for j = 1:ncrv
                    errSdnse_tru{j,i} = err_evl_tru( Sdnse_cell_i, j );
                    errSdnse_nse{j,i} = err_evl_nse( Sdnse_cell_i, j );
                    errSdnse_sys{j,i} = err_evl_sys( Sdnse_cell_i, j );

                    res_tru(:,i) = res_tru(:,i) + ...
                    sum( (  ...
                            Sdnse_cell_i{j}(2:end,:) ...
                            - pts_ref_cell{j}(2:end,:) ...
                    ).^2 ,2);
                    res_sys(i) = res_sys(i) + ...
                    sum((   ...
                        Sdnse_cell_i{j}((end-ndep+1):end,:) ...
                        - ode_.dnxu(Sdnse_cell_i{j}(1:(ndim-ndep),:)) ...
                        ).^2, 2);
                end
            end
            err_tru = sqrt(res_tru);
            [err_tru_min,i_err_tru_min] = min(err_tru,[],2)
            [iwrite_err_tru_min,err_rat_tru_min] = deal(iwrite_full(i_err_tru_min), err_tru_min./err_nse_ref)

            err_sys = sqrt(res_sys);
            [err_sys_min,i_err_sys_min] = min(err_sys,[],2)
            [iwrite_err_sys_min,err_rat_sys_min] = deal(iwrite_full(i_err_sys_min), err_sys_min./err_nse_sys)
            % iwrite_err_sys_min = iwrite_full(i_err_sys_min)

            [res_min,i_res_min] = min(dnse_summary.residuals);

            % icrv_plot = 1:6;
            icrv_plot = [1,3,6];

            plt_spc = plts_out(1);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',0.5,...
                                                        '.',12, ...
                                                        [0 0 0 1]);
            plt_spc = LD_plots.plot_pts(pts_ref_cell(icrv_plot), ...
                                    meta0,plt_spc,spc);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('none',0.5,...
                                                        '.',12, ...
                                                        [1 0 0 1]);
            plt_spc = LD_plots.plot_pts(pts_nse_cell(icrv_plot), ...
                                    meta0,plt_spc,spc);


            color_mat_short = flip(spring( nwrite_short ),1);

            color_mat_full = flip(cool( double(dnse_summary.nsmooth) ),1);
            % color_mat_long = flip(cool( nwrite_long ),1);
            color_mat_long = color_mat_full( iwrite_plot((nwrite_short+1):end) , :);

            color_mat = [color_mat_short ; color_mat_long];

            [spc.lspec,spc.lw,spc.mspec,spc.ms] = deal('-',1,...
                                                        'none',3);
            for i = 1:len_iwrite_plot
                spc.color = [color_mat(i,:) 1];
                plt_spc = LD_plots.plot_pts(Sdnse_cell(icrv_plot,i), ...
                                        meta0,plt_spc,spc);
            end
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                        'none',3, ...
                                                        [[color_mat(end,:)] 0.3] );
                                                        % [LD_plots.green5 1]  );
            plt_spc = LD_plots.plot_pts(Sobs_f(:,1), ...
                                    meta0,plt_spc,spc);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                        'none',3, ...
                                                        [[color_mat(end,:)] 1] );
                                                        % [LD_plots.green5 1]  );
            plt_spc = LD_plots.plot_pts(Sobs_f(icrv_plot,1), ...
                                    meta0,plt_spc,spc);


            plt_cnv = plts_out(2);
            [tdim1,tdim2] = deal(2,4);
            plt_cnv = plt_cnv.init_tiles_safe(tdim1,tdim2);
            hold(plt_cnv.axs, 'on');
            box(plt_cnv.axs,'on');
            axs_cnv = plt_cnv.axs;
            axs_mat_cnv = plt_cnv.axs_mat();

            plot_ax = @(ax_,x_,y_,s_) plot(ax_, ...
                x_, ...
                y_, ...
                'Color',s_.color, ...
                'Marker', s_.mspec, ...
                'LineStyle', s_.lspec, ...
                'MarkerSize',s_.ms, ...
                'LineWidth',s_.lw);

            axs_res = axs_mat_cnv(1,1);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                '.',8, ...
                                                                [1 0 0] );
            plot_ax(axs_res, ...
                1:(dnse_summary.nsmooth), ...
                dnse_summary.residuals, ...
                spc);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('none',1,...
                                                                '.',10, ...
                                                                [0 0 1]);
            plot_ax(axs_res, ...
                    dnse_summary.iwrite, ...
                    residuals_iwrite, ...
                    spc);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('none',1,...
                                                                's',6, ...
                                                                [0 0 1]);
            plot_ax(axs_res, ...
                    iwrite_plot, ...
                    residuals_iplot_iwrite, ...
                    spc);

            set(axs_res, ...
                'YScale','log', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_res, ...
                ['$$ i_{\mathrm{smooth}} $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_res, ...
                ['$$ \| \tilde{S}_{i} - \tilde{S}_{i-1} \| $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axs_err = axs_mat_cnv(1,2);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                '.',8, ...
                                                                [0 0 0] );
            plot_ax(axs_err, ...
                    iwrite_full, ...
                    err_tru(1,:), ... % /err_nse_ref(1)
                    spc);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                '.',8, ...
                                                                LD_plots.blue5 );
            plot_ax(axs_err, ...
                    iwrite_full, ...
                    err_tru(2,:), ... % /err_nse_ref(2)
                    spc);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                'o',3, ...
                                                                LD_plots.green4 );
            plot_ax(axs_err, ...
                    iwrite_full, ...
                    err_sys(1,:), ...% /err_nse_sys(1,:)
                    spc);

            set(axs_err, ...
                'YScale','log', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_err, ...
                ['$$ i_{\mathrm{smooth}} $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_err, ...
                ['$$ \| \tilde{S}_{i} - \hat{S} \| $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axs_rnk = axs_mat_cnv(1,3);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                '.',8, ...
                                                                [1 0 0] );
            plot_ax(axs_rnk, ...
                    1:(dnse_summary.nsmooth), ...
                    dnse_summary.ranks, ...
                    spc);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('none',1,...
                                                                '.',10, ...
                                                                [0 0 1]);
            plot_ax(axs_rnk, ...
                    dnse_summary.iwrite, ...
                    ranks_iwrite, ...
                    spc);

            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('none',1,...
                                                                's',6, ...
                                                                [0 0 1]);
            plot_ax(axs_rnk, ...
                    iwrite_plot, ...
                    ranks_iplot_iwrite, ...
                    spc);


            set(axs_rnk, ...
                'YScale','linear', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_rnk, ...
                ['$$ i_{\mathrm{smooth}} $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_rnk, ...
                ['$$ \mathrm{rank} ( \mathbf{R}_{i}^{(k)} ) $$'], 'Interpreter','Latex','FontSize',14);


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axs_sigma = axs_mat_cnv(2,3);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                '.',8, ...
                                                                [0 0 0] );
            plot_ax(axs_sigma, ...
                    1:ndof, ...
                    Rsvd_g_0_tru.s, ...
                    spc);

            spc.color = [1 0 0];
            plot_ax(axs_sigma, ...
                    1:ndof, ...
                    Rsvd_g_0.s, ...
                    spc);

            for i = 1:len_iwrite_plot
                spc.color = [color_mat(i,:) 1];
                plot(axs_sigma, ...
                    1:ndof, ...
                    Smat(:,i), ...
                    'Color',spc.color, ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end

            set(axs_sigma, ...
                'YScale','log', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_sigma, ...
                ['$$ i_{\theta} $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_sigma, ...
                ['$$ \sigma_{i_{\theta}} ( \mathbf{R}_{i}^{(k)} ) $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axs_theta = axs_mat_cnv(2,1);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                'none',1, ...
                                                                [0 0 0] );
            for i = 1:ndof
                plot(axs_theta, ...
                    1:ndof, ...
                    Rsvd_g_0_tru.V(:,i), ...
                    'Color',[spc.color, Rsvd_g_0_tru.s(end)/Rsvd_g_0_tru.s(i) ], ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end

            spc.color = [1 0 0];
            for i = 1:ndof
                plot(axs_theta, ...
                    1:ndof, ...
                    Rsvd_g_0.V(:,i), ...
                    'Color',[spc.color, Rsvd_g_0.s(end)/Rsvd_g_0.s(i) ], ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end

            spc.color = color_mat(1,:);
            for i = 1:ndof
                plot(axs_theta, ...
                    1:ndof, ...
                    Rsvd_g_s.V(:,i), ...
                    'Color',[spc.color, Rsvd_g_s.s(end)/Rsvd_g_s.s(i) ], ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end

            spc.color = color_mat(end,:);
            for i = 1:ndof
                plot(axs_theta, ...
                    1:ndof, ...
                    Rsvd_g_f.V(:,i), ...
                    'Color',[spc.color, Rsvd_g_f.s(end)/Rsvd_g_f.s(i) ], ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end

            set(axs_theta, ...
                'YScale','linear', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_theta, ...
                ['$$ i_{\theta} $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_theta, ...
                ['$$ \mathbf{v}_{i_{\theta}} ( \mathbf{R}_{i}^{(k)} ) $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            comp_spec_err = @(s_,V_) sqrt( sum( abs(SVT_g_0_tru*( s_ .* V_ )), 1));
            % comp_spec_err = @(s_,V_) sqrt( sum( (SVT_g_0_tru*( s_' .* V_ )).^2, 1));

            axs_ker = axs_mat_cnv(2,2);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                '.',8, ...
                                                                [0 0 0] );
            plot_ax(axs_ker, ...
                    1:ndof, ...
                    comp_spec_err(Rsvd_g_0_tru.s , Rsvd_g_0_tru.V), ...
                    spc);

            spc.color = [1 0 0];
            plot_ax(axs_ker, ...
                    1:ndof, ...
                    comp_spec_err(Rsvd_g_0.s , Rsvd_g_0.V), ...
                    spc);

            for i = 1:len_iwrite_plot
                spc.color = [color_mat(i,:) 1];
                plot(axs_ker, ...
                     1:ndof, ...
                     comp_spec_err(Smat(:,i)', Vtns(:,:,i)), ...
                    'Color',spc.color, ...
                    'Marker', spc.mspec, ...
                    'LineStyle', spc.lspec, ...
                    'MarkerSize',spc.ms,'LineWidth',spc.lw)
            end

            set(axs_ker, ...
                'YScale','log', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_ker, ...
                ['$$ i_{\theta} $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_ker, ...
                ['$$ \sigma_{i_{\theta}} ( \mathbf{R}_{i}^{(k)} ) $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axs_trj = axs_mat_cnv(1,4);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('-',1,...
                                                                'none',8, ...
                                                                [0 0 0] );
            for i = 1:length(icrv_plot)
                spc.color = [0 0 0 0.5];
                plot_ax(axs_trj, ...
                        pts_ref_cell{icrv_plot(i)}(1,:), ...
                        abs(errSdnse_tru{icrv_plot(i),1}(1,:) ), ...
                        spc);

                spc.color = [0 0 1 0.5];
                plot_ax(axs_trj, ...
                        pts_ref_cell{icrv_plot(i)}(1,:), ...
                        abs(errSdnse_tru{icrv_plot(i),1}(2,:) ), ...
                        spc);

                spc.color = [LD_plots.green4 0.5];
                plot_ax(axs_trj, ...
                        pts_ref_cell{icrv_plot(i)}(1,:), ...
                        abs(errSdnse_sys{icrv_plot(i),1}(1,:) ), ...
                        spc);
            end

            set(axs_trj, ...
                'YScale','log', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            xlabel(axs_trj, ...
                ['$$ x $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs_trj, ...
                ['$$ | \tilde{s}_i - s_i | $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            axs_ver = axs_mat_cnv(2,4);
            [spc.lspec,spc.lw,spc.mspec,spc.ms,spc.color] = deal('none',1,...
                                                                '.',8, ...
                                                                [0 0 0] );

            spc.color = [0 0 0 1];
            for i = 1:length(icrv_plot)
                plot_ax(axs_ver, ...
                        errSdnse_nse{icrv_plot(i),1}(2,:) , ...
                        errSdnse_tru{icrv_plot(i),1}(2,:) , ...
                        spc);
            end

            set(axs_ver, ...
                'YScale','linear', 'XScale','linear', ...
                'TickLabelInterpreter','Latex','FontSize',14);
            % xlabel(axs_ver, ...
            %     ['$$ x $$'], 'Interpreter','Latex','FontSize',14);
            % ylabel(axs_trj, ...
            %     ['$$ | \tilde{s}_i - s_i | $$'], 'Interpreter','Latex','FontSize',14);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            plts_out(1) = plt_spc;
            plts_out(2) = plt_cnv;

            dnse_data = dnse_summary;

            dnse_data.Sdnse_cell = Sdnse_cell;
        end
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

            % ndns_range = [ 1:5 , 100:100:900 , 999];

            % ndns_short_range = 1;
            % ndns_short_range = 1:2;
            % ndns_short_range = 1:5;

            % ndns_long_range = [];
            % ndns_long_range = (ndns_short_range(end)+1):1:10;
            % ndns_long_range = (ndns_short_range(end)+1):1:50;
            % ndns_long_range = (ndns_short_range(end)+1):1:3;
            % ndns_long_range = (6):2:50;
            % ndns_long_range = (10):10:500;
            % ndns_long_range = (100):100:500;
            % ndns_long_range = (200):200:1000;

            ndns_short_range = 1:5;
            % ndns_long_range = [];
            % ndns_long_range = (ndns_short_range(end)+1):1:50;
            % ndns_long_range = (ndns_short_range(end)+1):2:50;
            % ndns_long_range = 10:10:50;
            ndns_long_range = 100:100:1000;

            % ndns_short_range = 1;
            % ndns_short_range = 1:2;
            % ndns_long_range = [];

            ndns_range = [ndns_short_range ndns_long_range];

            len_ndns_range = length(ndns_range);

            pts_cell = cell([Sref_.ncrv,len_ndns_range]);

            color_mat = flip(cool(len_ndns_range),1);

            color_mat_short = flip(spring( length(ndns_short_range) ),1);
            color_mat_long = flip(cool( length(ndns_long_range) ),1);
            color_mat = [color_mat_short ; color_mat_long];

            for i = 1:len_ndns_range
                idns = ndns_range(i);
                pts_cell(:,i) = Snse_.read_Sobs_cell(['.jsol_Rk_' num2str(idns)]);
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
            % spc.color = [LD_plots.green5 1];
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
