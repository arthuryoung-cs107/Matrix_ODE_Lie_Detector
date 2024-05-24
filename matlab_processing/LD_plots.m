classdef LD_plots
    properties (Constant)
        screens = get(0,'MonitorPositions');

        grey1 = [178/250, 186/250, 187/250];
        grey2 = [131/250, 145/250, 146/250];
        grey3 = [97/250, 106/250, 107/250];
        grey4 = [66/250, 73/250, 73/250];
        grey5 = [20/100 20/100 20/100];

        purple1 = [102/250, 0/250, 102/250];
        purple2 = [153/250, 0/250, 153/250];
        purple3 = [204/250, 0/250, 204/250];
        purple4 = [250/250, 0/250, 250/250];
        purple5 = [250/250, 50/250, 250/250];

        orange1 = [255/255 90/255 0];
        orange2 = [255/255 123/255 0];
        orange3 = [255/255 165/255 0];
        orange4 = [255/255 208/255 0];
        orange5 = [255/255 229/255 0];

        green1 = [88/250, 214/250, 141/250];
        green2 = [40/250, 180/250, 99/250];
        green3 = [34/250, 153/250, 84/250];
        green4 = [25/250, 111/250, 61/250];
        green5 = [0, 1, 0];

        blue1 = [120/250, 150/250, 250/250];
        blue2 = [52/250, 152/250, 219/250];
        blue3 = [39/250, 97/250, 141/250];
        blue4 = [10/250, 50/250, 150/250];
        blue5 = [0, 0, 1];

        red1 = [236/250, 112/250, 99/250];
        red2 = [192/250, 57/250, 43/250];
        red3 = [146/250, 43/250, 33/250];
        red4 = [100/250, 30/250 , 22/250];
        red5 = [1, 0 , 0];

    end
    properties
        name;

        fig;
        axs;

        tile;
    end
    methods
        function obj = LD_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_)
            obj.name = name_;
            if (nargin == 2)
                posdim_specs_ = grid_dim_;
                obj = obj.set_screen_posdim(posdim_specs_);
            elseif (nargin == 4)
                obj = obj.set_screen_posdim(grid_dim_,tile_dim_,origin_tile_);
            elseif (nargin == 5)
                obj = obj.set_screen_posdim(grid_dim_,tile_dim_,origin_tile_,screen_);
            end
        end
        function obj_out = set_screen_posdim(obj,grid_dim_,tile_dim_,origin_tile_,screen_)
            if (nargin==5)
                grid_dim = grid_dim_;
                tile_dim = tile_dim_;
                origin_tile = origin_tile_;
                screen = screen_;
            elseif (nargin == 4)
                grid_dim = grid_dim_;
                tile_dim = tile_dim_;
                origin_tile = origin_tile_;
                screen = 1;
            elseif(nargin == 2)
                specs = grid_dim_;
                grid_dim = specs.grid_dim;
                tile_dim = specs.tile_dim;
                origin_tile = specs.origin_tile;
                screen = specs.screen;
            end

            if (screen>size(LD_plots.screens,1))
                screen_specs = LD_plots.screens(1,:);
            else
                screen_specs = LD_plots.screens(screen,:);
            end

            [o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1);
            dels_grid = (d_screen)./[grid_dim(2) grid_dim(1)];
            dels_tile = dels_grid.*[tile_dim(2) tile_dim(1)];

            oy_plot = floor(o_screen(2) + dels_grid(2)*(grid_dim(1) - origin_tile(1)));
            ox_plot = floor(o_screen(1) + dels_grid(1)*(origin_tile(2)-1));
            ly_plot = floor(dels_tile(2));
            lx_plot = floor(dels_tile(1));

            posdim_use = [ox_plot oy_plot lx_plot ly_plot];
            props_struct = LD_plots.make_posdim_plot_specs(obj.name,posdim_use);

            obj_out = obj;
            obj_out.fig = figure;
            for i=1:size(props_struct, 1)
                obj_out.fig.set(props_struct{i, 1}, props_struct{i, 2});
            end
        end
        function obj_out = init_tiles_safe(obj,tdim1_,tdim2_)
            if (length(obj.axs))
                obj_out = obj;
            else
                obj_out = obj;
                clf(obj_out.fig);
                obj_out.tile = tiledlayout(obj_out.fig,tdim1_,tdim2_);
                obj_out.tile.TileSpacing = 'compact';
                obj_out.tile.Padding = 'compact';

                tile_num = tdim1_*tdim2_;
                obj_out.axs = gobjects(tile_num, 1);
                for i=1:tile_num
                    obj_out.axs(i) = nexttile(obj_out.tile);
                end
            end
        end
    end
    methods (Static)
        function plt = plot_S_icrv_divergence(S_,svd_cmp_,icrv_,solspc_plot_,plt_)
            plt = plt_.init_tiles_safe(2,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            spc = LD_plots.make_default_plot_specs;
            spc.color = [1 0 0];
            spc.lw = 2;

            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,icrv_);

            pts_cell = S_.pts_cell;
            pts_mat_i = pts_cell{icrv_};
            pts_nrmlz = abs(pts_mat_i);
            pts_nrmlz(pts_mat_i==0) = 1.0;
            xvec = pts_mat_i(1,:);
            [npts,ncrv] = deal(size(pts_mat_i,2),S_.ncrv);
            [nrows,ncols,rho] = deal(svd_cmp_.nrow/ncrv,svd_cmp_.ncol,max(svd_cmp_.rvec));
            kappa = ncols - rho;
            nconstr_dim = nrows/npts;
            Ki = svd_cmp_.Vtns(:,(rho+1):end,icrv_);
            ATtns = reshape(svd_cmp_.matT,ncols,nrows,ncrv);

            crv_inds = 1:ncrv;
            ncrv_inds = crv_inds ~= icrv_;

            [rel_div_mat,inner_K_mat] = deal(nan(ncrv,npts));
            [rel_div_ics,rel_div_net,inner_K_net] = deal(nan(ncrv,1));
            for j = 1:ncrv
                rel_div_j = abs(pts_cell{j}-pts_mat_i)./pts_nrmlz;
                sum_rel_div_j = sum(rel_div_j,1);
                rel_div_mat(j,:) = sum_rel_div_j;
                rel_div_ics(j) = sum_rel_div_j(1);
                rel_div_net(j) = sum(sum_rel_div_j);

                inner_K_j_sumpts = reshape(sum(reshape(abs((ATtns(:,:,j))'*Ki),nconstr_dim,npts,kappa),1),npts,kappa);

                inner_K_mat(j,:) = sum(inner_K_j_sumpts,2);
                inner_K_net(j) = sum(inner_K_mat(j,:),2);
            end
            [~,i_sort_ics] = sort(rel_div_ics);
            [~,i_sort_net] = sort(rel_div_net);
            [~,i_sort_inn] = sort(inner_K_net);

            i_mindiv_ics = i_sort_ics(2);
            i_mindiv_net = i_sort_net(2);
            i_mininn_net = i_sort_inn(2);

            hard_red = LD_plots.red5;
            hard_green = LD_plots.green5;
            hard_blue = LD_plots.blue5;
            nice_green = LD_plots.green4;
            nice_purple = LD_plots.purple5;
            nice_orange = LD_plots.orange1;

            spc.color = hard_blue;
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_mindiv_ics);
            spc.color = hard_green;
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_mindiv_net);
            spc.color = nice_green;
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_mininn_net);

            mspc = 'none';
            ms = 1;
            lspc = '-';
            lw = 1;
            alpha = 0.2;

            plot(axs(1),rel_div_net(ncrv_inds),inner_K_net(ncrv_inds), ...
            'Marker','o','MarkerSize',2,'LineStyle','none','LineWidth',lw,'Color',hard_red);
            plot(axs(1),rel_div_net(i_mindiv_ics),inner_K_net(i_mindiv_ics), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',hard_blue);
            plot(axs(1),rel_div_net(i_mindiv_net),inner_K_net(i_mindiv_net), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',hard_green);
            plot(axs(1),rel_div_net(i_mininn_net),inner_K_net(i_mininn_net), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',nice_green);

            plot(axs(2),xvec,cumsum(rel_div_mat(ncrv_inds,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[hard_red alpha]);
            plot(axs(2),xvec,cumsum(rel_div_mat(i_mindiv_ics,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_blue);
            plot(axs(2),xvec,cumsum(rel_div_mat(i_mindiv_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_green);
            plot(axs(2),xvec,cumsum(rel_div_mat(i_mininn_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',nice_green);

            plot(axs(3),xvec,rel_div_mat(ncrv_inds,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[hard_red alpha]);
            plot(axs(3),xvec,rel_div_mat(i_mindiv_ics,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_blue);
            plot(axs(3),xvec,rel_div_mat(i_mindiv_net,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_green);
            plot(axs(3),xvec,rel_div_mat(i_mininn_net,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',nice_green);


            plot(axs(4),rel_div_ics(ncrv_inds),inner_K_net(ncrv_inds), ...
            'Marker','o','MarkerSize',2,'LineStyle','none','LineWidth',lw,'Color',hard_red);
            plot(axs(4),rel_div_ics(i_mindiv_ics),inner_K_net(i_mindiv_ics), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',hard_blue);
            plot(axs(4),rel_div_ics(i_mindiv_net),inner_K_net(i_mindiv_net), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',hard_green);
            plot(axs(4),rel_div_ics(i_mininn_net),inner_K_net(i_mininn_net), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',nice_green);


            plot(axs(5),xvec,cumsum(inner_K_mat(ncrv_inds,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[hard_red alpha]);
            plot(axs(5),xvec,cumsum(inner_K_mat(i_mindiv_ics,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_blue);
            plot(axs(5),xvec,cumsum(inner_K_mat(i_mindiv_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_green);
            plot(axs(5),xvec,cumsum(inner_K_mat(i_mininn_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',nice_green);

            plot(axs(6),xvec,inner_K_mat(ncrv_inds,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[hard_red alpha]);
            plot(axs(6),xvec,inner_K_mat(i_mindiv_ics,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_blue);
            plot(axs(6),xvec,inner_K_mat(i_mindiv_net,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_green);
            plot(axs(6),xvec,inner_K_mat(i_mininn_net,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',nice_green);

            xlabel(axs(1),['$$ \sum \mathrm{err} ( \mathbf{C}^j | \mathbf{C}^i ) $$'], 'Interpreter','Latex','FontSize',14);
            xlabel(axs(4),['$$ \mathrm{err} (s_0^j | s_0^i) $$'], 'Interpreter','Latex','FontSize',14);
            xlabel([axs(2:3) axs(5:6)],['$$ x $$'], 'Interpreter','Latex','FontSize',14);

            ylabel([axs(1) axs(4)],['$$ \sum | \mathbf{A}^j \mathbf{K}^i | $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(2),['$$ \int \mathrm{err} (s^j | s^i) dx $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(3),['$$ \mathrm{err} (s^j | s^i) $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(5),['$$ \int \sum | \mathbf{A}^j \mathbf{K}^i | dx $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(6),['$$ \sum | \mathbf{A}^j \mathbf{K}^i | $$'], 'Interpreter','Latex','FontSize',14);


            % set(axs(1),'YScale', 'log','XScale','log','TickLabelInterpreter','Latex','FontSize',12);
            set([axs(1) axs(4)],'YScale', 'linear','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);
            set([axs(2:3) axs(5:6)],'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);
        end
        function plt = plot_S_svds(Sarray_,plt_)
            plt = plt_.init_tiles_safe(2,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            S0 = Sarray_(1);
            nset = length(Sarray_);

            mspc = 'none';
            ms = 1;
            lspc = '-';
            lw = 1;
            colors = turbo(nset);

            mat_stats = @(mat_) deal(min(mat_,[],2),median(mat_,2),max(mat_,[],2), 1:size(mat_,1));

            [smin_cell,smed_cell,smax_cell,inds_cell] = deal(cell(nset,1));
            [sglb_cell,inmn_cell,inmd_cell,inmx_cell] = deal(cell(nset,1));
            for i = 1:nset
                % [smin_cell{i},smed_cell{i},smax_cell{i},inds_cell{i}] = mat_stats(Sarray_(i).smat);
                [smin_cell{i},smed_cell{i},smax_cell{i},inds_cell{i}] = mat_stats(Sarray_(i).smat ./ Sarray_(i).smat(1,:) );
                mat_i = Sarray_(i).matT';
                [~,si,Vi] = svd(mat_i,'econ','vector');
                % sglb_cell{i} = si;
                sglb_cell{i} = si/si(1);
                [inmn_cell{i},inmd_cell{i},inmx_cell{i},~] = mat_stats(abs(mat_i*Vi)');
            end


            for i = 1:nset
                leg1(i) = plot(axs(1),inds_cell{i},smin_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            for i = 1:nset
                leg2(i) = plot(axs(2),inds_cell{i},smed_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            for i = 1:nset
                leg3(i) = plot(axs(3),inds_cell{i},smax_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end

            for i = 1:nset
                leg4(i) = plot(axs(4),inds_cell{i},sglb_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            for i = 1:nset
                leg5(i) = plot(axs(5),inds_cell{i},inmd_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            for i = 1:nset
                leg6(i) = plot(axs(6),inds_cell{i},inmx_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end

            xlabel(axs(1:3),['$$ i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(1),['min $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(2),['med. $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(3),['max $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);

            ylabel(axs(4),['$$ \sigma^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(5),['med. $$ \mathbf{A}_g \mathbf{V}^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(6),['max $$ \mathbf{A}_g \mathbf{V}^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);

            legend(axs(5), leg1,'Location', 'SouthOutside', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));

            set(axs,'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);

        end
        function [plt,err_res] = plot_S_relative_errors(Sarray_,Sref_,plt_)
            plt = plt_.init_tiles_safe(1,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            nset = length(Sarray_);

            mspc = 'none';
            ms = 1;
            lspc_med = '-';
            lspc_max = ':';
            lw = 1;
            colors = turbo(nset);

            [ndim,nobs,nc,ndep] = deal(Sref_.ndim,Sref_.nobs,Sref_.ncrv,Sref_.ndep);
            np = nobs/nc;

            pts_mat_ref = Sref_.pts_mat;
            pts_tns_ref = reshape(pts_mat_ref,ndim,[],nc);
            x_vec = pts_tns_ref(1,:,1);

            pts_nrmlz = abs(pts_mat_ref);
            pts_nrmlz(pts_nrmlz==0) = 1.0;

            abs_rel_diff = nan(ndim,nobs,nset);
            for i = 1:nset
                abs_rel_diff_i = abs(pts_mat_ref-Sarray_(i).pts_mat)./(pts_nrmlz);
                abs_rel_diff(:,:,i) = abs_rel_diff_i;
            end

            re_net = reshape(sum(abs_rel_diff,1),np,nc,nset);
            re_u = reshape(sum(abs_rel_diff(2:(ndep+1),:,:),1),np,nc,nset);
            re_dnxu = reshape(sum(abs_rel_diff((end-ndep+1):end,:,:),1),np,nc,nset);

            xlabel(axs,['$$ x $$'], 'Interpreter','Latex','FontSize',14);
            set(axs,'XLim',[min(x_vec) max(x_vec)]);
            for i = 1:nset
                leg1(i) = plot(axs(1),x_vec,median(re_net(:,:,i),2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_med,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
                plot(axs(1),x_vec,max(re_net(:,:,i),[],2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_max,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            ylabel(axs(1),['med. $$ \mathrm{err} (\hat{s} | s ) $$'], 'Interpreter','Latex','FontSize',14);
            for i = 1:nset
                leg2(i) = plot(axs(2),x_vec,median(re_u(:,:,i),2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_med,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
                plot(axs(2),x_vec,max(re_u(:,:,i),[],2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_max,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            ylabel(axs(2),['med. $$ \mathrm{err} (\hat{u} | u ) $$'], 'Interpreter','Latex','FontSize',14);
            for i = 1:nset
                leg3(i) = plot(axs(3),x_vec,median(re_dnxu(:,:,i),2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_med,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
                plot(axs(3),x_vec,max(re_dnxu(:,:,i),[],2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_max,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            ylabel(axs(3),['med. $$ \mathrm{err} (d_{\hat{x}}^n \hat{u} | d_x^n u ) $$'], 'Interpreter','Latex','FontSize',14);

            set(axs(1:3),'YScale', 'log', 'XScale', 'linear', 'TickLabelInterpreter','Latex','FontSize',12);
            legend(axs(2), leg2,'Location', 'SouthOutside', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));

            err_res = struct(   'pts_nrmlz',pts_nrmlz, ...
                                'abs_rel_diff',abs_rel_diff);
        end
        function plt = plot_n1q1_solspc(pts_cell_,plt_,spc_)
            plt = plt_.init_tiles_safe(2,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');

            axs = plt.axs;

            dnames = {'x'; 'u'; 'd_x u'};
            dim_order = [   1,2; ...
                            1,3; ...
                            2,3];

            ncrv = length(pts_cell_);

            for i_plot = 1:size(dim_order,1)
                axi = axs(i_plot);
                dims_i = dim_order(i_plot,:);
                for i = 1:ncrv
                    plot(axi,pts_cell_{i}(dims_i(1),:),pts_cell_{i}(dims_i(2),:), ...
                    'Marker',spc_.mspec,'MarkerSize',spc_.ms,'LineStyle',spc_.lspec,'LineWidth',spc_.lw,'Color',spc_.color);
                end
                xlabel(axi,['$$' dnames{dims_i(1)} '$$'], 'Interpreter','Latex','FontSize',14);
                ylabel(axi,['$$' dnames{dims_i(2)} '$$'], 'Interpreter','Latex','FontSize',14);
            end
        end
        function plt = plot_n2q1_solspc(pts_cell_,plt_,spc_)
            plt = plt_.init_tiles_safe(2,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');

            axs = plt.axs;

            dnames = {'x'; 'u'; 'd_x u'; 'd^2_x u'};
            dim_order = [   1,2; ...
                            1,3; ...
                            1,4; ...
                            2,3; ...
                            2,4; ...
                            3,4];

            ncrv = length(pts_cell_);

            for i_plot = 1:size(dim_order,1)
                axi = axs(i_plot);
                dims_i = dim_order(i_plot,:);
                for i = 1:ncrv
                    plot(axi,pts_cell_{i}(dims_i(1),:),pts_cell_{i}(dims_i(2),:), ...
                    'Marker',spc_.mspec,'MarkerSize',spc_.ms,'LineStyle',spc_.lspec,'LineWidth',spc_.lw,'Color',spc_.color);
                end
                xlabel(axi,['$$' dnames{dims_i(1)} '$$'], 'Interpreter','Latex','FontSize',14);
                ylabel(axi,['$$' dnames{dims_i(2)} '$$'], 'Interpreter','Latex','FontSize',14);
            end
        end
        function specs_out = posdim_specs(grid_dim_,tile_dim_,origin_tile_,screen_)
            if (nargin==3)
                screen = 1;
            else
                screen = screen_;
            end
            if (screen>size(LD_plots.screens,1))
                screen_specs = LD_plots.screens(1,:);
            else
                screen_specs = LD_plots.screens(screen,:);
            end
            specs_out = struct( 'grid_dim',grid_dim_, ...
                                'tile_dim',tile_dim_, ...
                                'origin_tile',origin_tile_, ...
                                'screen',screen);
        end
        function specs_out = make_default_plot_specs()
            specs_out = struct( 'lspec', '-', ...
                                'mspec', 'none', ...
                                'ms', 1, ...
                                'lw', 0.5, ...
                                'color', [0 0 0]);
        end
        function struct_out = make_posdim_plot_specs(name_in_, pos_in_)
            struct_out = {'Name', name_in_; 'Renderer', 'painters'; 'Position', pos_in_;};
        end
        function plt = plot_solspc(S_,plt_,spc_,icrvs_)
            pts_cell_full = S_.pts_cell;
            if (nargin==2)
                spc = LD_plots.make_default_plot_specs;
                pts_cell_plot = pts_cell_full;
            elseif (nargin==3)
                if (isstruct(spc_))
                    spc = spc_;
                    pts_cell_plot = pts_cell_full;
                else
                    spc = LD_plots.make_default_plot_specs;
                    icrvs = spc_;
                    pts_cell_plot = pts_cell_full(icrvs);
                end
            else
                spc = spc_;
                pts_cell_plot = pts_cell_full(icrvs_);
            end
            switch S_.eor
                case 1
                    switch S_.ndep
                        case 1
                            plt = LD_plots.plot_n1q1_solspc(pts_cell_plot,plt_,spc);
                    end
                case 2
                    switch S_.ndep
                        case 1
                            plt = LD_plots.plot_n2q1_solspc(pts_cell_plot,plt_,spc);
                    end
            end
        end
    end
end

function label_out = fixlabel(label_in_)
    label_out = replace(label_in_,'_',' ');
end
