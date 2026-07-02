classdef apv_plots
    properties (Constant)

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

        view_mat = [45, 45; 1, 0; 0, 90; 90, 0 ; 45, 0; 70, 10; -20, 10; -220, 10];

    end
    properties
        name;

        fig;
        axs;

        tile;
    end
    methods
        function obj = apv_plots(name_,grid_dim_,tile_dim_,origin_tile_,screen_)
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
        function axs_mat_out = axs_mat(obj)
            [tdim1,tdim2] = deal(obj.tile.GridSize(1),obj.tile.GridSize(2));
            axs_mat_out = (reshape(obj.axs,tdim2,tdim1))';
        end
        function obj_out = set_screen_posdim(obj,grid_dim_,tile_dim_,origin_tile_,screen_)
            if (nargin==2)
                % specs = grid_dim_;
                % grid_dim = specs.grid_dim;
                % tile_dim = specs.tile_dim;
                % origin_tile = specs.origin_tile;
                % screen = specs.screen;

                posdim_use = grid_dim_;
            else
                if (nargin==5)
                    grid_dim = grid_dim_;
                    tile_dim = tile_dim_;
                    origin_tile = origin_tile_;
                    screen = screen_;
                else (nargin == 4)
                    grid_dim = grid_dim_;
                    tile_dim = tile_dim_;
                    origin_tile = origin_tile_;
                    screen = 1;
                end

                sys_screens = apv_plots.get_sys_screens();
                % sys_screens = get(groot,'MonitorPositions');

                if (screen>size(sys_screens,1))
                    screen_specs = sys_screens(1,:); % default to screen 1
                else
                    screen_specs = sys_screens(screen,:);
                end

                [o_screen d_screen] = deal(screen_specs(1:2), screen_specs(3:4)-1);
                dels_grid = (d_screen)./[grid_dim(2) grid_dim(1)];
                dels_tile = dels_grid.*[tile_dim(2) tile_dim(1)];

                oy_plot = floor(o_screen(2) + dels_grid(2)*(grid_dim(1) - origin_tile(1)));
                ox_plot = floor(o_screen(1) + dels_grid(1)*(origin_tile(2)-1));
                ly_plot = floor(dels_tile(2));
                lx_plot = floor(dels_tile(1));

                posdim_use = [ox_plot oy_plot lx_plot ly_plot];
            end

            props_struct = apv_plots.make_posdim_plot_specs(obj.name,posdim_use);

            obj_out = obj;
            obj_out.fig = figure('MenuBar', 'none', 'ToolBar', 'none');
            for i=1:size(props_struct, 1)
                % obj_out.fig.set(props_struct{i, 1}, props_struct{i, 2});
                set(obj_out.fig,props_struct{i, 1},props_struct{i, 2})
            end
        end
        function obj_out = init_tiles_safe(obj,tdim1_,tdim2_)
            if (nargin==2)
                tdims = tdim1_;
                tdim1_=tdims(1);
                tdim2_=tdims(2);
            end
            if (~isempty(obj.axs))
                obj_out = obj;
            else
                obj_out = obj;
                clf(obj_out.fig);
                % figure(obj_out.fig);
                % obj_out.tile = subplot(tdim1_,tdim2_);
                obj_out.tile = tiledlayout(obj_out.fig,tdim1_,tdim2_,'TileSpacing','compact','Padding','compact');

                tile_num = tdim1_*tdim2_;
                obj_out.axs = gobjects(tile_num, 1);
                for i=1:tile_num
                    obj_out.axs(i) = nexttile(obj_out.tile);
                end
            end
        end
        function obj_out = set_subplots_safe(obj,dims_)
            if (~isempty(obj.axs))
                obj_out = obj;
            else
                obj_out = obj;
                clf(obj_out.fig);

                nax = size(dims_,1);
                obj_out.axs = gobjects(nax,1);

                for i = 1:nax
                    obj_out.axs(i) = subplot(dims_(i,1),dims_(i,2),dims_(i,3));
                end
            end
        end
        function obj_out = set_axes_safe(obj,dims_)
            if (~isempty(obj.axs))
                obj_out = obj;
            else
                obj_out = obj;
                clf(obj_out.fig);

                nax = size(dims_,1);
                obj_out.axs = gobjects(nax,1);

                for i = 1:nax
                    obj_out.axs(i) = axes('Position',dims_(i,:));
                end
            end
        end
        function show_toolbar(obj)
            set(obj.fig, 'ToolBar', 'Figure');
        end
        function show_menubar(obj)
            set(obj.fig, 'MenuBar', 'Figure');
        end
        function hide_toolbar(obj)
            set(obj.fig, 'ToolBar', 'none');
        end
        function hide_menubar(obj)
            set(obj.fig, 'MenuBar', 'none');
        end
        function axis_lims_out = get_axis_lims(obj)
            naxes = length(obj.axs);
            axis_lims_out = cell(naxes,1);
            for i = 1:naxes
                axis_lims_out{i} = axis(obj.axs(i));
            end
        end
        function set_axis_lims(obj,lims_)
            naxes = length(obj.axs);
            for i = 1:naxes
                axis(obj.axs(i),lims_{i});
            end
        end

        function filename = write_figure(obj,type_,dir_,name_)
            if (nargin == 4)
                name = name_;
            else
                name = obj.name;
            end

            filename = [dir_ name '.' type_];

            if (strcmp(type_,'pdf'))
                exportgraphics(obj.fig,filename,'ContentType','vector');
            else
                exportgraphics(obj.fig,filename);
            end
        end
    end

    methods (Static)
        function [plt,plt_1D] = plot_LDsol_model_summary(p_,mod_,d_)
            plt_jspc = apv_plots('jspc_model_summary', ...
                [2 1],...
                [1 1],[1 1], ...
                1);
            plt_jspc.show_toolbar

            ndep = length(mod_.fspace_0.Omap_b(:))-1;
            ndim = size(mod_.Smat,1);
            kor = ((ndim-1)/ndep) - 1;
            [plt_jspc,pspc] = apv_plots.plot_Sobs(plt_jspc,mod_.Sobs,d_);
            plt = plt_jspc;
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            tau_uN_RN1_tns = reshape(mod_.tau_uN_RN1_net(:,1,:),ndep,kor,[]);
            for i = 1:ndep
                tau_uiN = reshape(tau_uN_RN1_tns(i,:,:),kor,[]);
                for k = 2:(kor+1)
                    plot(axs_mat(i,k), ...
                    mod_.Smat(1,:), tau_uiN(k-1,:), ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', pspc.Color, ...
                    'Marker', pspc.Marker, ...
                    'MarkerSize', ceil(0.25*pspc.MarkerSize), ...
                    'MarkerFaceColor', [0 1 0], ...
                    'MarkerEdgeColor', [0 0 0] ...
                    );
                end
            end

            plt = p_;
            plt = plt.init_tiles_safe(1,2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            set(axs, ...
                'XScale', 'linear',  ...
                'TickLabelInterpreter','Latex', ...
                'FontSize',16 );
            xlabel(axs, '$$ i $$', 'Interpreter','Latex','FontSize',16);
            % plt.show_menubar
            plt.show_toolbar

            function plot_labelled_svds(axi_,svds_,lbls_,leg_loc_)
                nsvd = length(svds_(:));
                cmat = hsv(nsvd);
                if (nsvd==1)
                    leg_i(1) = plot(axi_, ...
                    1:length(svds_.s), svds_.s, ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', cmat(1,:), ...
                    'Marker', 's', ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', cmat(1,:), ...
                    'MarkerEdgeColor', cmat(1,:), ...
                    'DisplayName', lbls_{1} ...
                    );
                    % 'MarkerEdgeColor', [0 0 0], ...
                else
                    for i = 1:nsvd
                        leg_i(i) = plot(axi_, ...
                        1:length(svds_{i}.s), svds_{i}.s, ...
                        'LineStyle', 'none', ...
                        'LineWidth', 0.5, ...
                        'Color', cmat(i,:), ...
                        'Marker', 's', ...
                        'MarkerSize', 6, ...
                        'MarkerFaceColor', cmat(i,:), ...
                        'MarkerEdgeColor', cmat(i,:), ...
                        'DisplayName', lbls_{i} ...
                        );
                        % 'MarkerEdgeColor', [0 0 0], ...
                    end
                end
                set(axi_, ...
                    'YScale', 'log' );
                ylabel(axi_, '$ \sigma_i $', 'Interpreter','Latex','FontSize',16);
                legend(axi_, leg_i(1:nsvd),'Location', leg_loc_, 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',14);
            end
            rlbl = @(s_) [ ', \rho=' num2str(s_.r) '/' num2str(s_.dim) ];
            nlbl = @(s_) [ ', \| \cdot \|_* =' num2str(sum(s_.s/s_.s(1)),'%.1f') '/' num2str(s_.dim) ];

            % title(axi, '$$ \textrm{SVDs over parameter space, } \mathbb{R}^C $$', ...
            axi = axs(1);
            title(axi, ['SVDs over parameter space '], ...
                'Interpreter','Latex','FontSize',12 );
            svds_i = { ...
                mod_.H_N1_svd , '$H^N'; ...
                mod_.DprN_svd, '$D^N'; ...
                mod_.Rsvd_N1 , '$R^N'; ...
                mod_.Rsvd_N1_net, '$R^N_{\mathrm{net}}'; ...
                mod_.Tsvd_N1,  '$T^N'; ...
                mod_.Gsvd_N1,  '$G^N'; ...
                mod_.Gsvd_N1_net,  '$G^N_{\mathrm{net}}'; ...
                mod_.Gsvd_N1_TV, '$G^N \textrm{ (non-trivial) }'; ...
                mod_.Gsvd_N1_TV_net, '$G^N_{\mathrm{net}} \textrm{ (non-trivial) }'; ...
            };
            svds = cell([ size(svds_i,1),1 ]);
            labels = cell([length(svds),1]);
            for isvd = 1:size(svds_i,1)
                labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) nlbl(svds_i{isvd,1}) '$'];
                svds{isvd} = svds_i{isvd,1};
            end
            plot_labelled_svds(axi,svds,labels,'EastOutside');

            % title(axi, '$$ \textrm{SVDs over jet space, } \mathbb{R}^B $$', ...
            axi = axs(2);
            title(axi, ['SVDs over jet space '], ...
                'Interpreter','Latex','FontSize',12 );
            svds_i = { ...
                mod_.GN1_sO_basis , '$\Lambda^{(0)} |_{s_O} W_G '; ...
                mod_.GN1_sNO_basis, '$\Lambda^{(N)} |_{s_O} W_G'; ...
            };
            svds = cell([ size(svds_i,1),1 ]);
            labels = cell([length(svds),1]);
            for isvd = 1:size(svds_i,1)
                labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) nlbl(svds_i{isvd,1}) '$'];
                svds{isvd} = svds_i{isvd,1};
            end
            plot_labelled_svds(axi,svds,labels,'EastOutside');

        end
        %% apv plots
        function [plt,pspc] = plot_1D_model_summary(p1_,p2_,mod_,d_)
            plt = p1_;
            plt_svds = p2_;

            % plt.show_menubar
            plt.show_toolbar
            plt_svds.show_toolbar

            mods_1D = mod_.jspc_1Dmods;
            ndep = length( mods_1D(:) );
            kor_obs = size(mods_1D{1}.Smat,1) - 2;
            ncrv = length( mod_.Sobs(:) );
            Smat = mod_.Smat;
            xvec = Smat(1,:);
            utns = reshape(Smat(2:end,:),ndep,kor_obs+1,[]);
            nobs = length(xvec);

            cmat_crv = hsv(ncrv);
            cmat_ord = cool(kor_obs);

            % nplt_adtl = 0;
            % plt = plt.init_tiles_safe(ndep,1+kor_obs+nplt_adtl);
            [plt,pspc] = apv_plots.plot_Sobs(plt,mod_.Sobs,d_);
            % [plt,pspc] = apv_plots.plot_Sobs(plt,mod_.Sobs{1},d_);
            axs_mat = plt.axs_mat;
            lims0 = plt.get_axis_lims();

            plt_svds = plt_svds.init_tiles_safe(ndep,3);
            hold(plt_svds.axs, 'on');
            box(plt_svds.axs,'on');
            axs_svds_mat = plt_svds.axs_mat;

            Plen_1D_crv = mods_1D{1}.fspace.Plen;
            ntheta_1D_crv = 2*Plen_1D_crv;
            for i = 1:ndep
                mod_i = mods_1D{i};
                bor_i = mod_i.fspace.bor;
                crv_svds_i = mod_i.RNsvds_crv;
                crv_sub_svds_i = mod_i.RN_crv_sub_svds;
                glb_svds_i = mod_i.Rksvds_glb;
                RN_i_name = ['R^{(N)}_{' num2str(i) '}'];
                crv_c_name = @(c_) ['\textrm{ (crv' num2str(c_) '), }'];

                rank_c_name = @(c_) ['\rho=' num2str(crv_svds_i{c_}.r) '/' num2str(crv_svds_i{c_}.dim)];
                rank_bc_name = @(b_,c_) ['\rho=' num2str(crv_sub_svds_i{b_,c_}.r) '/' num2str(crv_sub_svds_i{b_,c_}.dim)];

                svd_name_crv_i = @(c_) ['$$' crv_c_name(c_) rank_c_name(c_) '$$'];
                svd_name_sub_crv_i = @(b_,c_) ['$$' crv_c_name(c_) 'b=' num2str(b_) ', ' rank_bc_name(b_,c_) '$$'];

                Rk_i_name = @(k_) ['R^{' num2str(k_) '}_{' num2str(i) '}'];
                rank_k_name = @(k_) ['\rho=' num2str(glb_svds_i{k_}.r) '/' num2str(glb_svds_i{k_}.dim)];
                svd_name_glb_k = @(k_) ['$$' Rk_i_name(k_) ',' rank_k_name(k_) ' \textrm{ (glb) } $$'];

                tau_k_i = @(k_) reshape(mod_i.tau_uN_crv(k_,:),size(xvec));
                % bsub = min([10,bor_i])
                % bsub = 3
                bsub = bor_i;
                tau_k_i_alt = @(k_) reshape(mod_i.tau_uN_crv_sub(k_,:,bsub),size(xvec));
                % tau_k_i_alt = @(k_) reshape(mod_i.tau_uk_glb(k_,:),size(xvec));

                % xvec, reshape(mod_i.tau_uN_crv(k-1,:),size(xvec)), ...
                for k = 2:(kor_obs+1)

                    if (ndep == 1)
                        plot(axs_mat(i,k), ...
                        xvec, tau_k_i(k-1), ...
                        'LineStyle', 'none', ...
                        'LineWidth', 0.5, ...
                        'Color', pspc.Color, ...
                        'Marker', pspc.Marker, ...
                        'MarkerSize', ceil(0.75*pspc.MarkerSize), ...
                        'MarkerFaceColor', [1 0 0], ...
                        'MarkerEdgeColor', [0 0 0] ...
                        );
                        plot(axs_mat(i,k), ...
                        xvec, tau_k_i_alt(k-1), ...
                        'LineStyle', 'none', ...
                        'LineWidth', 0.5, ...
                        'Color', pspc.Color, ...
                        'Marker', pspc.Marker, ...
                        'MarkerSize', ceil(0.5*pspc.MarkerSize), ...
                        'MarkerFaceColor', [0 0 1], ...
                        'MarkerEdgeColor', [0 0 0] ...
                        );
                    end

                    leg_glb_k(k-1) = plot(axs_svds_mat(i,1), ...
                    1:length(mod_i.Rksvds_glb{k-1}.s), mod_i.Rksvds_glb{k-1}.s, ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Marker', 's', ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', cmat_ord(k-1,:), ...
                    'MarkerEdgeColor', [0 0 0], ...
                    'DisplayName', svd_name_glb_k(k-1) ...
                    );
                end
                legend(axs_svds_mat(i,1), leg_glb_k(1:kor_obs),'Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',12);
                % plot(axs_mat(i,end), ...
                % xvec, mod_i.tau_uk_glb(end,:), ...
                % 'LineStyle', 'none', ...
                % 'LineWidth', 0.5, ...
                % 'Color', pspc.Color, ...
                % 'Marker', pspc.Marker, ...
                % 'MarkerSize', ceil(0.5*pspc.MarkerSize), ...
                % 'MarkerFaceColor', apv_plots.orange1, ...
                % 'MarkerEdgeColor', [0 0 0] ...
                % );

                for icrv = 1:ncrv
                    leg_crv_i(icrv) = plot(axs_svds_mat(i,2), ...
                    1:ntheta_1D_crv, mod_i.RNsvds_crv{icrv}.s, ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Marker', 's', ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', cmat_crv(icrv,:), ...
                    'MarkerEdgeColor', [0 0 0], ...
                    'DisplayName', svd_name_crv_i(icrv) ...
                    );
                    for b = 1:bor_i
                        leg_sub_crv_i(b,icrv) = plot(axs_svds_mat(i,3), ...
                        1:length(mod_i.RN_crv_sub_svds{b,icrv}.s), mod_i.RN_crv_sub_svds{b,icrv}.s, ...
                        'LineStyle', '-', ...
                        'LineWidth', 0.5, ...
                        'Color', cmat_crv(icrv,:), ...
                        'Marker', 'none', ...
                        'MarkerSize', pspc.MarkerSize, ...
                        'MarkerFaceColor', cmat_crv(icrv,:), ...
                        'MarkerEdgeColor', [0 0 0], ...
                        'DisplayName', svd_name_sub_crv_i(b,icrv) ...
                        );
                    end
                end
                title(axs_svds_mat(i,2), ['$$' RN_i_name '$$'], ...
                    'Interpreter','Latex','FontSize',12 );
                title(axs_svds_mat(i,3), ['$$' RN_i_name ' \textrm{ (sub, } b=' num2str(bsub) ') $$'], ...
                    'Interpreter','Latex','FontSize',12 );
                legend(axs_svds_mat(i,2), leg_crv_i(1:ncrv),'Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',12);
                legend(axs_svds_mat(i,3), leg_sub_crv_i(bsub,:),'Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',12);
            end
            set(axs_svds_mat(:), ...
                'XScale', 'linear',  ...
                'YScale', 'log',  ...
                'TickLabelInterpreter','Latex', ...
                'FontSize',10 );
            ylabel(axs_svds_mat(:), '$$ \sigma_i $$', 'Interpreter','Latex','FontSize',16);
            xlabel(axs_svds_mat(:), '$$ i $$', 'Interpreter','Latex','FontSize',16);

            plt.set_axis_lims(lims0);

        end
        function [plt,plt_1D] = plot_model_summary(p_,mod_,d_)
            plt_jspc = apv_plots('jspc_model_summary_1D', ...
                [2 1],...
                [1 1],[1 1], ...
                1);
            plt_1D_svds = apv_plots('svds_model_summary_1D', ...
                [2 1],...
                [1 1],[2 1], ...
                1);
            [plt_1D,pspc] = apv_plots.plot_1D_model_summary(plt_jspc,plt_1D_svds,mod_,d_);
            plt = plt_1D;
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            tau_uN_tns = reshape(mod_.jspc_N1mod.tau_uN_net(:,1,:),size(axs_mat,1),size(axs_mat,2)-1,[]);
            for i = 1:size(axs_mat,1)
                % tau_uiN = reshape(mod_.tau_uN_crv(i,:,:),size(axs_mat,2),[]);
                tau_uiN = reshape(tau_uN_tns(i,:,:),size(axs_mat,2)-1,[]);
                for k = 2:size(axs_mat,2)
                    % mod_.Smat(1,1:343), tau_uiN(k-1,1:343), ...
                    plot(axs_mat(i,k), ...
                    mod_.Smat(1,:), tau_uiN(k-1,:), ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', pspc.Color, ...
                    'Marker', pspc.Marker, ...
                    'MarkerSize', ceil(0.25*pspc.MarkerSize), ...
                    'MarkerFaceColor', [0 1 0], ...
                    'MarkerEdgeColor', [0 0 0] ...
                    );
                end
            end

            plt = p_;
            plt = plt.init_tiles_safe(2,2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            set(axs, ...
                'XScale', 'linear',  ...
                'TickLabelInterpreter','Latex', ...
                'FontSize',16 );
            xlabel(axs, '$$ i $$', 'Interpreter','Latex','FontSize',16);
            % plt.show_menubar
            plt.show_toolbar

    axi = axs(1);
            % svds_ = mvp_jspc_model.vertcat_glb_SVDs(mod_);
            Rqksvd_glb = mod_.Rqksvd_glb_cell;
            Rksvd_glb = mod_.Rksvds_glb;
            Rksvd_sub_glb = mod_.Rk_glb_sub_svds;

            [ndep,kor_obs] = size(Rqksvd_glb);
            bor = size(Rksvd_sub_glb,1);

            % svds_ = vertcat( Rqksvd_glb(:),Rksvd_glb(:),Rksvd_sub_glb(:) );
            svds = vertcat( Rqksvd_glb(:),Rksvd_glb(:) );
            klbl = @(s_) [ ', \kappa=' num2str( s_.dim-s_.r ) '/' num2str(s_.dim) ];
            labels = cell([length(svds),1]);
            isvd = 1;
            for k = 1:kor_obs
                for i = 1:ndep
                    labels{isvd} = ['$R^{' num2str(k) '}_{' num2str(i) '}' klbl(svds{isvd}) '\textrm{ (glb) } $'] ;
                    isvd = isvd+1;
                end
            end
            for k = 1:kor_obs
                labels{isvd} = ['$R^{' num2str(k) '}' klbl(svds{isvd}) '\textrm{ (glb) } $'] ;
                isvd = isvd+1;
            end
            % for k = 1:kor_obs
            %     for b = 1:bor
            %     % for b = [1, bor]
            %         labels{isvd} = ['$R^{' num2str(k) '}' klbl(svds{isvd}) '\textrm{ (glb,b=' num2str(b) ') $'];
            %         isvd = isvd+1;
            %     end
            % end
            leg_loc = 'EastOutside'; % leg_loc = 'NorthEast';
            % leg_ncols = min([length(svds(:)),2]);
            leg_ncols = 1;
            function plot_labelled_svds(svds_)
                nsvd = length(svds_(:));
                cmat = hsv(nsvd);
                if (nsvd==1)
                    leg_i(1) = plot(axi, ...
                    1:length(svds_.s), svds_.s, ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', cmat(1,:), ...
                    'Marker', 's', ...
                    'MarkerSize', 6, ...
                    'MarkerFaceColor', cmat(1,:), ...
                    'MarkerEdgeColor', cmat(1,:), ...
                    'DisplayName', labels{1} ...
                    );
                    % 'MarkerEdgeColor', [0 0 0], ...
                else
                    for i = 1:nsvd
                        leg_i(i) = plot(axi, ...
                        1:length(svds_{i}.s), svds_{i}.s, ...
                        'LineStyle', 'none', ...
                        'LineWidth', 0.5, ...
                        'Color', cmat(i,:), ...
                        'Marker', 's', ...
                        'MarkerSize', 6, ...
                        'MarkerFaceColor', cmat(i,:), ...
                        'MarkerEdgeColor', cmat(i,:), ...
                        'DisplayName', labels{i} ...
                        );
                        % 'MarkerEdgeColor', [0 0 0], ...
                    end
                end
                set(axi, ...
                    'YScale', 'log' );
                ylabel(axi, '$$ \sigma_i $$', 'Interpreter','Latex','FontSize',16);
                legend(axi, leg_i(1:nsvd),'Location', leg_loc, 'Interpreter', 'Latex', 'NumColumns',leg_ncols,'FontSize',14);
            end
            plot_labelled_svds(svds);

            Rqksvd_crv = mod_.Rqksvds_crv_cell; % mod_.Rqksvds_crv_cell(:,end,:)
            RNsvds_crv = mod_.RNsvds_crv;
            RNsvds_sub_crv = mod_.RN_crv_sub_svds; % mod_.RN_crv_sub_svds(end,:)
            ncrv = length(RNsvds_crv);
            % svds_ = vertcat( Rqksvd_crv(:),RNsvds_crv(:),RNsvds_sub_crv(:) );

    axi = axs(2);
            svds = RNsvds_crv(:);
            labels = cell([length(svds),1]);

            klbl = @(s_) [ ', \kappa=' num2str( s_.dim-s_.r ) '/' num2str(s_.dim) ];
            isvd = 1;
            for j = 1:ncrv
                labels{isvd} = ['$$R^{(N)}' klbl(svds{isvd}) '\textrm{ (crv=' num2str(j) ') } $$'] ;
                isvd = isvd+1;
            end
            leg_loc = 'EastOutside'; % leg_loc = 'EastOutside';
            % leg_ncols = min([length(svds(:)),2]);
            leg_ncols = 1;
            plot_labelled_svds(svds);

            vthx_N_glb = mod_.vthx_N_glb_svd;
            lvs_glb = mod_.lvs_svd;
            LamN_glb = mod_.LamN_svd;
            Gsvd = mod_.Gsvd;
            Hsvd_0 = mod_.Hsvd_0;
            Hsvd_Nm1 = mod_.jspc_N1mod.Hsvd_Nm1;

            Lam_dkxuq_glb = mod_.Lam_dkxuq_svd_glb_cell;
            Lam_dkxu_glb = mod_.Lam_dkxu_svd_glb;
            RN_YLdNm1xu_glb = mod_.RN_YLdNm1xu_glb_svd;

    axi = axs(3);
            svds_i = { ...
                vthx_N_glb , '$\vartheta^x |_{R^N}'; ...
                lvs_glb , '$\lambda |_{\{(x,u)_j\}}'; ...
                LamN_glb , '$\Lambda^{(N)} |_{\{(x,u^{(N)})_j\}}'; ...
                Gsvd , '$G'; ...
                Hsvd_0 , '$H^0'; ...
                Hsvd_Nm1 , '$H^{N-1}'; ...
            };
            rlbl = @(s_) [ ', \rho=' num2str(s_.r) '/' num2str(s_.dim) ];
            nlbl = @(s_) [ ', \| \cdot \|_* =' num2str(sum(s_.s/s_.s(1)),'%.1f') '/' num2str(s_.dim) ];

            svds = cell([ size(svds_i,1),1 ]);
            labels = cell([length(svds),1]);
            for isvd = 1:size(svds_i,1)
                labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) nlbl(svds_i{isvd,1}) '$'];
                svds{isvd} = svds_i{isvd,1};
            end
%             isvd = size(svds_i,1)+1;
%             for k = 1:kor_obs
%                 for i = 1:ndep
%                     svds_{isvd} = Lam_dkxuq_glb{i,k};
% prefix_i = ['\Lambda_{d_x^{' num2str(k) '} u_{' num2str(i) '}} |_{(x,u^{(' num2str(k) ')})_j}'];
% labels{isvd} = ['$$' prefix_i ', \rho = ' num2str(svds_{isvd}.r)  '$$'];
%                     isvd = isvd + 1;
%                 end
%             end
            leg_loc = 'EastOutside'; % leg_loc = 'EastOutside';
            leg_ncols = 1;
            plot_labelled_svds(svds);
            set(axi, ...
                'YScale', 'log' ); %'YScale', 'linear' );

            RN1_glb_svd = mod_.jspc_N1mod.Rsvd_glb;
            DprN_svd = mod_.jspc_N1mod.DprN_svd;
            RN1_net_svd = mod_.jspc_N1mod.Rsvd_net;
            vth_RN1_net_svd = mod_.jspc_N1mod.vth_net_svd;
            GNm1_glb_svd = mod_.jspc_N1mod.GNm1_glb_svd;
            GNm1_net_svd = mod_.jspc_N1mod.GNm1_net_svd;

    axi = axs(4);
            svds_i = { ...
                RN1_glb_svd , '$R^{(N)}_{\textrm{glb}}' ; ...
                DprN_svd , '$D^{(N)}_{\textrm{glb}}' ; ...
                RN1_net_svd , '$R^{(N)}_{\textrm{net}}' ; ...
                vth_RN1_net_svd , '$\vartheta_{R^{(N)}_{\textrm{net}}}' ; ...
                GNm1_glb_svd , '$G^{(1)}_{\textrm{glb}}' ; ...
                GNm1_net_svd , '$G^{(1)}_{\textrm{net}}' ; ...
            };
            rlbl = @(s_) [ ', \rho=' num2str(s_.r) '/' num2str(s_.dim) ];
            nlbl = @(s_) [ ', \| \cdot \|_* =' num2str(sum(s_.s/s_.s(1)),'%.1f') '/' num2str(s_.dim) ];

            svds = cell([ size(svds_i,1),1 ]);
            labels = cell([length(svds),1]);
            for isvd = 1:size(svds_i,1)
                labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) nlbl(svds_i{isvd,1}) '$'];
                svds{isvd} = svds_i{isvd,1};
            end
            leg_loc = 'EastOutside'; % leg_loc = 'EastOutside';
            leg_ncols = 1;
            plot_labelled_svds(svds);

            % RN_crv_Yvthxglb = mod_.RNsvds_crv_Yvthxglb;
            % cmat = cool(vthx_N_glb.r);
            % % for i = flip(1:(vthx_N_glb.r))
            % for i = 1:(vthx_N_glb.r)
            %     leg_i(i) = plot(axi, ...
            %     1:length(vthx_N_glb.s), vthx_N_glb.V(:,i), ...
            %     'LineStyle', '-', ...
            %     'LineWidth', 0.5, ...
            %     'Color', cmat(i,:), ...
            %     'Marker', 'none', ...
            %     'MarkerSize', 6, ...
            %     'MarkerFaceColor', cmat(i,:), ...
            %     'MarkerEdgeColor', [0 0 0], ...
            %     'DisplayName', ['$$ \mathbf{v}_{' num2str(i) '}$$'] ...
            %     );
            %     % fspc.print_vshort_polynomial_theta_z(vthx_N_glb.V(:,i),mod_.fdat.Pmat,'x');
            % end
            % legend(axi, leg_i(1:vthx_N_glb.r),'Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns',leg_ncols,'FontSize',10);
            % ylabel(axi, '$$ \mathbf v_{i,j} $$', 'Interpreter','Latex','FontSize',16);

            % apv_plots.set_containing_axis_lims(axs(:));

            % pause
        end
        function plt = plot_SVDs(p_,svds_,d_)
            plt = p_;

            plt = plt.init_tiles_safe(1,2);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            set(axs, ...
                'XScale', 'linear',  ...
                'TickLabelInterpreter','Latex', ...
                'FontSize',16 );
            xlabel(axs, '$$ i $$', 'Interpreter','Latex','FontSize',16);

            nsvd = length(svds_(:));
            cmat = hsv(nsvd);

            for i = 1:nsvd
                leg1(i) = plot(axs(1), ...
                1:length(svds_{i}.s), svds_{i}.s, ...
                'LineStyle', 'none', ...
                'LineWidth', 0.5, ...
                'Color', cmat(i,:), ...
                'Marker', 's', ...
                'MarkerSize', 6, ...
                'MarkerFaceColor', cmat(i,:), ...
                'MarkerEdgeColor', [0 0 0], ...
                'DisplayName', ['SVD' num2str(i)] ...
                );
            end
            set(axs(1), ...
                'YScale', 'log' );
            ylabel(axs(1), '$$ \sigma_i $$', 'Interpreter','Latex','FontSize',16);

            legend(axs(1), leg1,'Location', 'EastOutside', 'Interpreter', 'Latex');
            % legend(axs(1), leg1,'Location', 'SouthWest', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));

            for i = 1:nsvd
                plot(axs(2), ...
                1:length(svds_{i}.s), svds_{i}.V', ...
                'LineStyle', '-', ...
                'LineWidth', 0.5, ...
                'Color', cmat(i,:), ...
                'Marker', 'none' ...
                );
                % 'Marker', 's', ...
                % 'MarkerSize', 6, ...
                % 'MarkerFaceColor', cmat(i,:), ...
                % 'MarkerEdgeColor', [0 0 0] ...
            end
            % alpha(axs(2),0.5)
            ylabel(axs(2), '$$ \mathbf{v}_{j,i} $$', 'Interpreter','Latex','FontSize',16);

        end
        function [plt,pspc] = plot_Sobs(p_,S_,d_)
            plt = p_;

            pspc = apv_plots.verify_plotspecs(d_);
            pspc.MarkerSize = 8;
            % pspc.MarkerSize = 4;

            eor = d_.eor;
            ndep = d_.ndep;
            ndim = 1+ndep*(eor+1);

            if (isempty(plt.axs))
                if (ndim==3)
                    plt = plt.init_tiles_safe(ndep,eor+2);
                else
                    plt = plt.init_tiles_safe(ndep,eor+1);
                end
            end
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            axs_set = axs_mat(:,1:(eor+1));
            set(axs_set, ...
                'YScale', 'linear',  ...
                'XScale', 'linear',  ...
                'TickLabelInterpreter','Latex', ...
                'FontSize',16 );
            xlabel(axs_set, '$$ x $$', 'Interpreter','Latex','FontSize',16);

            if (iscell(S_))
                Smat = ldaux.Scell_2_Smat(S_,ndim);

            elseif ( length(size(S_)) == 2 )
                Smat = S_;

            else
                Smat = reshape( S_, ndim, [] );
            end

            xvec = Smat(1,:);
            utns = reshape(Smat(2:end,:),ndep,eor+1,[]);

            yname = @(i_,k_) ['$$ d_x^' num2str(k_) 'u_' num2str(i_)  '$$'];

            for k = 1:(eor+1)
                for i = 1:ndep
                    plot(axs_mat(i,k), ...
                    xvec, reshape(utns(i,k,:),size(xvec)), ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', pspc.Color, ...
                    'Marker', pspc.Marker, ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', pspc.Color, ...
                    'MarkerEdgeColor', [0 0 0] ...
                    );
                    ylabel(axs_mat(i,k), yname(i,k-1), ...
                     'Interpreter','Latex','FontSize',16);
                end
            end

            %% if eor=ndep=1, plot 3D jetspace surface
            if (ndim==3)
                plot3( axs(3), xvec, Smat(2,:), Smat(3,:) , ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', pspc.Color, ...
                    'Marker', pspc.Marker, ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', pspc.Color, ...
                    'MarkerEdgeColor', [0 0 0] ...
                    );
                zlabel(axs(3), '$$ d^1_x u_1 $$', 'Interpreter','Latex','FontSize',16);
                view(axs(3), apv_plots.view_mat(6, :));
            end

            if ( isfield(d_,'LineStyle') )
                if ( ~strcmp( d_.LineStyle , 'none' ) )

                    if (iscell(S_))
                        ncrv = length(S_);
                        Scell = S_;
                    elseif ( length(size(S_)) == 2 )
                        ncrv = 1;
                        Scell = cell([1,1]);
                        Scell{1} = S_;
                    else
                        ncrv = size(S_,3);
                        Scell = cell([ncrv,1]);
                        for j = 1:ncrv
                            Scell{j} = S_(:,:,j);
                        end
                    end
                    iv = 2;
                    for k = 1:(eor+1)
                        for i = 1:ndep
                            for j = 1:ncrv
                                plot(axs_mat(i,k), ...
                                Scell{j}(1,:), Scell{j}(iv,:), ...
                                'LineStyle', d_.LineStyle, ...
                                'LineWidth', 0.5, ...
                                'Color', pspc.Color, ...
                                'Marker', 'none' ...
                                );
                            end
                            iv = iv+1;
                        end
                    end
                    if (ndim==3)
                        for j = 1:ncrv
                            plot3( axs(3), ...
                            Scell{j}(1,:), Scell{j}(2,:), Scell{j}(3,:), ...
                            'LineStyle', d_.LineStyle, ...
                            'LineWidth', 0.5, ...
                            'Color', pspc.Color, ...
                            'Marker', 'none' ...
                            );
                        end
                    end
                end
            end
        end

        %% true plotting utilities

        function spc_out = verify_plotspecs(spc_)
            spc_out = spc_;

            if (~isfield(spc_out,'Color'))
                spc_out.Color = [0 0 0];
            end
            if (~isfield(spc_out,'Marker'))
                spc_out.Marker = 'o';
            end
            if (~isfield(spc_out,'LineStyle'))
                spc_out.Color = [0 0 0];
            end

        end

        function dim_out = near_squaredim(num_)

            root_floor0 = floor(sqrt(double(num_)));
            root_floor = root_floor0;

            while ( rem(num_,root_floor) ~= 0 )

                root_floor = root_floor - 1;

                if (root_floor==1)
                    break;
                end
            end

            if (root_floor==1)
                dim_out = [root_floor0+1,root_floor0+1];
            else
                dim_alt = num_/root_floor;
                dims = [root_floor dim_alt];
                dim_out = [min(dims) max(dims)];
            end
        end
        function set_containing_cscale_lims(axs_)
            axs = axs_(:);
            naxes = length(axs);
            cslim_mat = nan(naxes,2);
            for i = 1:naxes
                cslim_mat(i,:) = clim(axs(i));
            end

            clim_set = [min(cslim_mat(:,1)),max(cslim_mat(:,2))];

            for i = 1:naxes
                clim(axs(i),clim_set);
            end

        end
        function set_containing_axis_lims(axs_)
            axs = axs_(:);
            naxes = length(axs);
            axlim_tns = nan( 2, size( axis(axs(1)),2)/2, naxes );
            for i = 1:naxes
                axlim_tns(:,:,i) = reshape(axis(axs(i)),2,[]);
            end

            min_axlw = reshape(min(axlim_tns(1,:,:),[],3),1,[]);
            max_axhi = reshape(max(axlim_tns(2,:,:),[],3),1,[]);

            axlim_cnt = reshape([min_axlw;max_axhi],1,[]);

            for i = 1:naxes
                axis(axs(i), axlim_cnt);
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

        %% meta

        function sys_screens_out = get_sys_screens()
            sys_screens_out = get(groot,'MonitorPositions');
            if (~ismac) % works fine if osx
                arch = getenv('ARCH');
                [istart,iend] = regexp(arch,'mac');
                if ( (length(istart)*length(iend)) == 0 ) % works fine if osx
                    [istart,iend] = regexp(arch,'win'); % assume works fine if windows
                    if ( (length(istart)*length(iend)) == 0 ) % assume linux
                        [~,raw_out] = system('xrandr -q'); % query xrandr for linux displays
                        [match,nomatch] = regexp(raw_out,'\w*connected\w*','match','split');
                        nmonitors = length(match);

                        screen0 = nomatch{1,1};
                        i_c = regexp(screen0,'current');
                        i_m = regexp(screen0,'maximum');
                        numstrings0 = extract(screen0(i_c:(i_m-1)),digitsPattern);
                        height0 = str2num(numstrings0{2,1});

                        sys_screens_out = nan(nmonitors,4);
                        for i = 1:nmonitors
                            substr_i = nomatch{1,i+1};
                            i_p = regexp(substr_i,'(');
                            dimstring = substr_i(2:(i_p-1));
                            numstrings = extract(dimstring,digitsPattern);
                            heighti_true = str2num(numstrings{2});

                            if (length(regexp(dimstring,'primary')))
                                heighti = floor(0.9*heighti_true);
                            else
                                heighti = heighti_true;
                            end

                            sys_screens_out(i,:) = [    str2num(numstrings{3}), ...
                                                        height0-str2num(numstrings{4})-heighti, ...
                                                        str2num(numstrings{1}), ...
                                                        heighti];
                        end
                    end
                end
            end
        end
    end
end
