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
        function obj = LD_plots(name_)
            obj.name = name_;
        end
        function obj_out = set_screen_posdim(obj,grid_dim_,tile_dim_,origin_tile_,screen_)
            if (nargin==4)
                screen = 1;
            else
                screen = screen_;
            end

            if (screen>size(LD_plots.screens,1))
                screen_specs = LD_plots.screens(1,:);
            else
                screen_specs = LD_plots.screens(screen,:);
            end

            [o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1);
            dels_grid = (d_screen)./[grid_dim_(2) grid_dim_(1)];
            dels_tile = dels_grid.*[tile_dim_(2) tile_dim_(1)];

            oy_plot = floor(o_screen(2) + dels_grid(2)*(grid_dim_(1) - origin_tile_(1)));
            ox_plot = floor(o_screen(1) + dels_grid(1)*(origin_tile_(2)-1));
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
        function plt = plot_n1q1_solspc(S_,plt_)
            if (nargin<3)
                plt = LD_plots([S_.dat_name '_n1q1_solspc']);
                plt = plt.set_screen_posdim([4 4], [1 4], [4 1], 1);
            else
                plt = plt_;
            end
            plt = plt.init_tiles_safe(2,3);
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
            if (nargin<3)
                plt = LD_plots([S_.dat_name '_n2q1_solspc']);
                plt = plt.set_screen_posdim([4 4], [1 4], [1 1], 1);
            else
                plt = plt_;
            end
            plt = plt.init_tiles_safe(2,3);
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
    end
end
