classdef LD_plots
    properties (Constant)
        screens = get(0,'MonitorPositions');
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
        function plt = plot_S_svds(Sarray_,plt_)
            plt = plt_.init_tiles_safe(1,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            S0 = Sarray_(1);
            nset = length(Sarray_);

            mspc = 'none';
            ms = 1;
            lspc_min = '-';
            lspc_med = '-';
            lspc_max = '-';
            lw = 1;
            colors = turbo(nset);

            mat_stats = @(mat_) deal(min(mat_,[],2),median(mat_,2),max(mat_,[],2), 1:size(mat_,1));

            xlabel(axs,['$$ i $$'], 'Interpreter','Latex','FontSize',14);
            for i = 1:nset
                leg1(i) = plot(axs(1),1:(Sarray_(i).ndof),min(Sarray_(i).smat,[],2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_min,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            ylabel(axs(1),['min $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);

            for i = 1:nset
                leg2(i) = plot(axs(2),1:(Sarray_(i).ndof),median(Sarray_(i).smat,2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_med,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            ylabel(axs(2),['med. $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);

            for i = 1:nset
                leg3(i) = plot(axs(3),1:(Sarray_(i).ndof),max(Sarray_(i).smat,[],2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc_max,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(Sarray_(i).dat_name));
            end
            ylabel(axs(3),['max $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);

            legend(axs(2), leg2,'Location', 'SouthOutside', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));

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
        function plt = plot_n1q1_solspc(S_,plt_)
            plt = plt_.init_tiles_safe(2,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');

            axs = plt.axs;

            dnames = {'x'; 'u'; 'd_x u'};
            dim_order = [   1,2; ...
                            1,3; ...
                            2,3];

            spc = LD_plots.make_default_plot_specs;

            ncrv = S_.ncrv;
            pts_cell = S_.pts_cell;

            for i_plot = 1:size(dim_order,1)
                axi = axs(i_plot);
                dims_i = dim_order(i_plot,:);
                for i = 1:ncrv
                    plot(axi,pts_cell{i}(dims_i(1),:),pts_cell{i}(dims_i(2),:), ...
                    'Marker',spc.mspec,'MarkerSize',spc.ms,'LineStyle',spc.lspec,'LineWidth',spc.lw,'Color',spc.color);
                end
                xlabel(axi,['$$' dnames{dims_i(1)} '$$'], 'Interpreter','Latex','FontSize',14);
                ylabel(axi,['$$' dnames{dims_i(2)} '$$'], 'Interpreter','Latex','FontSize',14);
            end
        end
        function plt = plot_n2q1_solspc(S_,plt_)
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

            spc = LD_plots.make_default_plot_specs;

            ncrv = S_.ncrv;
            pts_cell = S_.pts_cell;

            for i_plot = 1:size(dim_order,1)
                axi = axs(i_plot);
                dims_i = dim_order(i_plot,:);
                for i = 1:ncrv
                    plot(axi,pts_cell{i}(dims_i(1),:),pts_cell{i}(dims_i(2),:), ...
                    'Marker',spc.mspec,'MarkerSize',spc.ms,'LineStyle',spc.lspec,'LineWidth',spc.lw,'Color',spc.color);
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
        function plt = plot_solspc(S_,plt_)
            switch S_.eor
                case 1
                    switch S_.ndep
                        case 1
                            LD_plots.plot_n1q1_solspc(S_,plt_);
                    end
                case 2
                    switch S_.ndep
                        case 1
                            LD_plots.plot_n2q1_solspc(S_,plt_);
                    end
            end
        end
    end
end

function label_out = fixlabel(label_in_)
    label_out = replace(label_in_,'_',' ');
end
