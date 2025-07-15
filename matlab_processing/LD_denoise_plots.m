classdef LD_denoise_plots < LD_plots
    properties

    end
    methods
        function obj = LD_denoise_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_)
            obj@LD_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_);
        end
    end
    methods (Static)
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
            set(axs(1:3),'YScale','log', 'XScale','linear');
            % set(axs(3),'YScale','log', 'XScale','log');

            % LD_plots.set_containing_axis_lims(axs_mat(1,1))
            % LD_plots.set_containing_axis_lims(axs_mat(2,1))

        end
        function plt = plot_curve_estimates(Splt_,plt_,spc_,Sref_,Snse_,icrv_)
            plt = plt_;
            [tdim1,tdim2] = deal(3,3);
            plt = plt.init_tiles_safe(tdim1,tdim2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = (reshape(axs,tdim2,tdim1))';

            [eor,ndep] = deal(Sref_.eor,Sref_.ndep);
            [crvs_ref_,crvs_nse_] = deal(Sref_.crvs,Snse_.crvs);

            i_crv = icrv_;
            kor_u_hat = eor;

            crv_i = crvs_ref_(i_crv);
            jet_i = crv_i.jets;
            jet_in = crv_i.make_jets(kor_u_hat);
            % x_check = 0:0.5:25;
            % x_check = jet_i.pts_mat(1,:);
            % x_check = jet_i.xh_vec;
            x_check = sort([jet_i.xh_vec,jet_i.pts_mat(1,:)]);
            inds_01 = 1:(2):length(x_check);
            inds_h = 2:(2):(length(x_check)-1);

            % u_check = jet_i.u_hat(x_check)
            [u_check,dxu_check] = jet_i.u_hat(x_check,kor_u_hat);
            [u_checkn,dxu_checkn] = jet_in.u_hat(x_check,kor_u_hat);

            x_check_n = x_check;

            sigma_unp1_ref = ones(eor+2,ndep);
            lam_dnp2xu_ref = ones(ndep,1);

            sigma_unp1 = [reshape(Snse_.sigma_s(2:end),ndep,[])' ; reshape(Snse_.sigma_dnp1xu,1,[])];
            % lam_dnp2xu = 1e-3*ones(ndep,1);
            lam_dnp2xu = ones(ndep,1);

            crv_i_n = crvs_nse_(i_crv);

            jet_i_n = crv_i_n.jets;
            rjet_i_n = LD_jets.regularize_jet(jet_i_n,sigma_unp1_ref,lam_dnp2xu_ref);
            rsjet_i_n = LD_jets.regularize_jet(jet_i_n,sigma_unp1,lam_dnp2xu);
            jet_in_n = crv_i_n.make_jets(kor_u_hat);

            [u_check_n,dxu_check_n] = jet_i_n.u_hat(x_check_n,eor);
            [u_rcheck_n,dxu_rcheck_n] = rjet_i_n.u_hat(x_check_n,eor);
            [u_rscheck_n,dxu_rscheck_n] = rsjet_i_n.u_hat(x_check_n,eor);
            [u_checkn_n,dxu_checkn_n] = jet_in_n.u_hat(x_check_n,eor);

            absres = @(u_,y_) abs(u_-y_);
            abserr = @(u_,y_) abs((u_-y_)./u_);
            mdplt = @(x_) ones(1,2)*median(x_);
            medabserr = @(u_,y_) abs((u_-y_)./u_);
            row_stats = @(row_) deal(min(row_,[],2),median(row_,2),mean(row_,2),max(row_,[],2));

            meta0 = LD_observations_set.make_meta_data(crv_i.eor,crv_i.ndep);
            spc = spc_;
            spc.color = LD_plots.green5;
            LD_plots.plot_pts([x_check; u_check; dxu_check],meta0,Splt_,spc);
            spc.color = LD_plots.orange1;
            LD_plots.plot_pts([x_check_n; u_check_n; dxu_check_n],meta0,Splt_,spc);
            % spc.color = LD_plots.orange4;
            % LD_plots.plot_pts([x_check_n; u_checkn_n; dxu_checkn_n],meta0,Splt_,spc);
            spc.color = LD_plots.green4;
            LD_plots.plot_pts([x_check_n; u_rcheck_n; dxu_rcheck_n],meta0,Splt_,spc);
            spc.color = LD_plots.blue3;
            LD_plots.plot_pts([x_check_n; u_rcheck_n; dxu_rcheck_n],meta0,Splt_,spc);

            [ucll,dxucll,d2xucll] = deal(cell([3,1]));
            [ucll{1},ucll{2},ucll{3}] = deal(   u_check_n, ...
                                                u_rcheck_n, ...
                                                u_rscheck_n);
            [dxucll{1},dxucll{2},dxucll{3}] = deal( dxu_check_n(1,:), ...
                                                    dxu_rcheck_n(1,:), ...
                                                    dxu_rscheck_n(1,:));
            [d2xucll{1},d2xucll{2},d2xucll{3}] = deal(  dxu_check_n(2,:), ...
                                                        dxu_rcheck_n(2,:), ...
                                                        dxu_rscheck_n(2,:));
            ctns = ones(2,3,3);
            ctns(:,:,1) = [LD_plots.orange1;LD_plots.orange4];
            ctns(:,:,2) = [LD_plots.green4;LD_plots.green1];
            ctns(:,:,3) = [LD_plots.blue3;LD_plots.blue1];

            icll = cell([2,1]);
            icll{1} = inds_01;
            icll{2} = inds_h;

            evl = {abserr ; absres};

            for iax = 1:2
                axi = axs_mat(iax,1);
                for i = 1:size(ctns,3)
                    for j = 1:length(icll)
                        plot(axi, ...
                            x_check(icll{j}),evl{iax}(u_check(icll{j}),ucll{i}(icll{j})), ...
                            's','Color',ctns(j,:,i),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                        plot(axi, ...
                            [x_check(1),x_check(end)],mdplt(evl{iax}(u_check(icll{j}),ucll{i}(icll{j}))), ...
                            '-','Color',ctns(j,:,i),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                    end
                end
            end
            for iax = 1:2
                axi = axs_mat(iax,2);
                for i = 1:size(ctns,3)
                    for j = 1:length(icll)
                        plot(axi, ...
                            x_check(icll{j}),evl{iax}(dxu_check(1,icll{j}),dxucll{i}(icll{j})), ...
                            's','Color',ctns(j,:,i),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                        plot(axi, ...
                            [x_check(1),x_check(end)],mdplt(evl{iax}(dxu_check(1,icll{j}),dxucll{i}(icll{j}))), ...
                            '-','Color',ctns(j,:,i),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                    end
                end
            end
            for iax = 1:2
                axi = axs_mat(iax,3);
                for i = 1:size(ctns,3)
                    for j = 1:length(icll)
                        plot(axi, ...
                            x_check(icll{j}),evl{iax}(dxu_check(2,icll{j}),d2xucll{i}(icll{j})), ...
                            's','Color',ctns(j,:,i),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                        plot(axi, ...
                            [x_check(1),x_check(end)],mdplt(evl{iax}(dxu_check(2,icll{j}),d2xucll{i}(icll{j}))), ...
                            '-','Color',ctns(j,:,i),'MarkerSize',spc.ms,'LineWidth',spc.lw)
                    end
                end
            end

            axi = axs_mat(3,1);
            plot(axi,  ...
                abs(u_check_n(inds_01)-u_check(inds_01)), ...
                abs(u_check_n(inds_01)-u_rcheck_n(inds_01)), ...
                's','Color',LD_plots.green4,'MarkerSize',spc.ms,'LineWidth',spc.lw)
            plot(axi,  ...
                abs(u_check_n(inds_01)-u_check(inds_01)), ...
                abs(u_check_n(inds_01)-u_rscheck_n(inds_01)), ...
                's','Color',LD_plots.blue3,'MarkerSize',spc.ms,'LineWidth',spc.lw)

            axi = axs_mat(3,2);
            plot(axi,  ...
                abs(dxu_check_n(1,inds_01)-dxu_check(1,inds_01)), ...
                abs(dxu_check_n(1,inds_01)-dxu_rcheck_n(1,inds_01)), ...
                's','Color',LD_plots.green4,'MarkerSize',spc.ms,'LineWidth',spc.lw)
            plot(axi,  ...
                abs(dxu_check_n(1,inds_01)-dxu_check(1,inds_01)), ...
                abs(dxu_check_n(1,inds_01)-dxu_rscheck_n(1,inds_01)), ...
                's','Color',LD_plots.blue3,'MarkerSize',spc.ms,'LineWidth',spc.lw)

            axi = axs_mat(3,3);
            plot(axi,  ...
                abs(dxu_check_n(2,inds_01)-dxu_check(2,inds_01)), ...
                abs(dxu_check_n(2,inds_01)-dxu_rcheck_n(2,inds_01)), ...
                's','Color',LD_plots.green4,'MarkerSize',spc.ms,'LineWidth',spc.lw)
            plot(axi,  ...
                abs(dxu_check_n(2,inds_01)-dxu_check(2,inds_01)), ...
                abs(dxu_check_n(2,inds_01)-dxu_rscheck_n(2,inds_01)), ...
                's','Color',LD_plots.blue3,'MarkerSize',spc.ms,'LineWidth',spc.lw)

            % [min_res01,med_res01,avg_res01,max_res01] = row_stats(absres(u_check(inds_01),u_check_n(inds_01)));
            % [min_res_h,med_res_h,avg_res_h,max_res_h] = row_stats(absres(u_check(inds_h),u_check_n(inds_h)));

            axs_set = axs_mat(:);
            % set(axs_mat(1,:), ...
            set(axs_set, 'TickLabelInterpreter','Latex','FontSize',12);
            % set(axs_mat(:),'YScale','linear', 'XScale','linear');
            set(axs_mat(1:2,:),'YScale','log', 'XScale','linear');
            set(axs_mat(3,:),'YScale','linear', 'XScale','linear');
            % set(axs_mat(3,:),'YScale','log', 'XScale','log');

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
            spc.color = [0 0 1 1];
            plt_out = LD_plots.plot_solspc(Sref_,plt_out,spc,icrv_plot);

            [spc.mspec,spc.ms] = deal('s',4);
            [spc.lspec,spc.lw] = deal('none',1);
            spc.color = [1 0 0 1];
            plt_out = LD_plots.plot_solspc(Snse_,plt_out,spc,icrv_plot);
        end
    end
end
