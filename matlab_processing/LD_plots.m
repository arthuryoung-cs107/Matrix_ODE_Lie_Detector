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

        view_mat = [45, 45; 1, 0; 0, 90; 90, 0 ; 45, 0; 70, 10; -20, 10; -220, 10];

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
            set(obj_out.fig, 'MenuBar', 'none', 'ToolBar', 'none');
            for i=1:size(props_struct, 1)
                obj_out.fig.set(props_struct{i, 1}, props_struct{i, 2});
            end
        end
        function obj_out = init_tiles_safe(obj,tdim1_,tdim2_)
            if (~isempty(obj.axs))
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
    end
    methods (Static)
        function plt = plot_global_component(S_,svd_cmp_,svd_cmp2_,svd_Lmat_,icrv_in_,fspc_,solspc_plot_,plt_)
            plt = plt_.init_tiles_safe(1,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            hard_red = LD_plots.red5;
            hard_green = LD_plots.green5;
            hard_blue = LD_plots.blue5;
            nice_green = LD_plots.green4;
            nice_purple = LD_plots.purple1;
            nice_pink = LD_plots.purple5;
            nice_orange = LD_plots.orange1;

            blue_mat = LD_plots.get_blue_mat;
            green_mat = LD_plots.get_green_mat;
            orange_mat = LD_plots.get_orange_mat;

            spc = LD_plots.make_default_plot_specs;
            spc.lw = 2;

            Amat_full = svd_cmp_.matT';
            [ncol,ncrv] = deal(svd_cmp_.ncol,length(svd_cmp_.rvec));
            Atns_full = permute(reshape(svd_cmp_.matT,ncol,[],ncrv),[2 1 3]);
            [s_glb,V_glb] = deal(svd_cmp_.svec_glb,svd_cmp_.Vmat_glb);

            [eor,ndep] = deal(S_.eor,S_.ndep);
            meta_ = S_.meta_data;
            nvar = ndep+1;
            pts_cell = S_.pts_cell;

            % % ncmp = 50;
            % % ncmp = 10;
            % % ncmp = 9;
            % ncmp = 21;
            % % ncmp = ncol;
            % K_glb = V_glb(:,(end-ncmp+1):end);
            % % K_glb = V_glb(:,(end-ncmp+1));
            % absAK = abs(pagemtimes(Atns_full,K_glb));
            % % absAK_evl = reshape(sum(sum(absAK,1),2),ncrv,1);
            % % absAK_evl = reshape(max(sum(absAK,2),[],1),ncrv,1);
            % absAK_evl = reshape(median(sum(absAK,2),1),ncrv,1);
            % % absAK_evl = reshape(max(max(absAK,[],1),[],2),ncrv,1);
            % [~,isort_absAK_evl] = sort(absAK_evl);

            Atns_cmp = Atns_full;
            % Atns_cmp = permute(reshape(svd_cmp_.matT ./ S_.pts_mat(end,:) ,svd_cmp_.ncol,[],ncrv),[2 1 3]);
            % Atns_cmp = permute(pagemtimes(svd_cmp_.Vtns, ...
            %                     eye(svd_cmp_.ncol).*reshape(svd_cmp_.smat,svd_cmp_.ncol,1,ncrv)), ...
            %                     [2 1 3]);
            % Atns_cmp = permute(pagemtimes(svd_cmp_.Vtns, ...
            %                     eye(svd_cmp_.ncol).*reshape(svd_cmp_.smat ./ svd_cmp_.smat(1,:), ...
            %                     svd_cmp_.ncol,1,ncrv)),[2 1 3]);

            % Atns_cmp = permute(reshape(svd_cmp_.matT ./ sqrt(sum(svd_cmp_.matT.*svd_cmp_.matT,1)), ...
            %                         ncol,[],ncrv),[2 1 3]);

            % Atns_cmp = permute(reshape(svd_cmp2_.matT,svd_cmp2_.ncol,[],ncrv),[2 1 3]);
            % Atns_cmp = permute(pagemtimes(svd_cmp2_.Vtns, ...
            %                     eye(svd_cmp2_.ncol).*reshape(svd_cmp2_.smat,svd_cmp2_.ncol,1,ncrv)), ...
            %                     [2 1 3]);
            % Atns_cmp = permute(reshape(pagemtimes(permute(svd_cmp2_.Vtns,[2 1 3]), ...
            %                                 eye(svd_cmp2_.ncol).*reshape(svd_cmp2_.smat,svd_cmp2_.ncol,1,ncrv)), ...
            %                                 ncol,[],ncrv),[2 1 3]);

            ncol_A = size(Atns_cmp,2);
            Vtns_cmp = svd_cmp2_.Vtns;
            if (size(Vtns_cmp,1)<ncol_A)
                Vtns_cmp = LD_aux.Ytns_Vtns_mult(svd_cmp2_.Ytns_L,svd_cmp2_.Vtns);
            end

            icrv_ = icrv_in_;
            inds_crv = 1:ncrv;
            inds_not_icrv = inds_crv(inds_crv ~= icrv_);

            rvec_cmp = svd_cmp2_.rvec;
            [min_rho,max_rho,rho_icrv] = deal(min(rvec_cmp),max(rvec_cmp),rvec_cmp(icrv_));
            ncol_cmp = size(Vtns_cmp,2);

            % % rho_use = max_rho
            % % rho_use = min_rho;
            % rho_use = rho_icrv;
            % ik_use = (rho_use+1):ncol_cmp;

            % ndim_use = 5;
            % ndim_use = (ncol_cmp-max_rho)-1;
            % ndim_use = floor((ncol_cmp-max_rho)/5);
            ndim_use = floor((ncol_cmp-rho_icrv)/2);

            % Amat_not_icrv = reshape(permute(Atns_cmp(:,:,inds_not_icrv),[2 1 3]),ncol_A,[])';
            Amat_not_icrv = Amat_full;

            % [~,s_glb_rstr,V_glb_rstr] = svd(Amat_not_icrv*Vtns_cmp(:,(rho_icrv+1):end,icrv_),'econ','vector');
            % k_glb_rstr = V_glb_rstr(:,end);
            % [~,isort_mag_k_glb_rstr] = sort(abs(k_glb_rstr));
            % ik_use = (isort_mag_k_glb_rstr((end-ndim_use+1):end)+double(rho_icrv));
            % AV_raw = pagemtimes(Atns_cmp,Vtns_cmp);
            % % AK_raw = AV_raw(:,ik_use,:);
            % K_glb = Vtns_cmp(:,ik_use,icrv_);

            rho_use = max_rho;
            [~,~,V_glb_not_icrv] = svd(Amat_not_icrv,'econ','vector');
            AK_glb = Vtns_cmp(:,(rho_use+1):end,icrv_)'*V_glb_not_icrv(:,(rho_use+1):end);
            AK_rstr_glb = Vtns_cmp(:,(rho_use+1):end,icrv_)*AK_glb;
            AV_raw = pagemtimes(Atns_cmp,AK_rstr_glb);
            K_glb = AK_rstr_glb;

            fro_AK_Vglb = norm(AK_rstr_glb,'fro')
            fro_Vglb = norm(V_glb_not_icrv(:,(rho_use+1):end),'fro')
            fro_AK_Vglb_diff = norm(AK_rstr_glb-V_glb_not_icrv(:,(rho_use+1):end),'fro')
            fro_AK_Vglb_diff/fro_Vglb

            row_stats = @(row_) deal(min(row_,[],2),median(row_,2),mean(row_,2),max(row_,[],2));
            AK_rstr_glb_colmags = sqrt(sum(AK_rstr_glb.*AK_rstr_glb,1));
            [min_colmag,med_colmag,avg_colmag,max_colmag] = row_stats(AK_rstr_glb_colmags)

            [~,imn_colmag] = min(AK_rstr_glb_colmags,[],2)
            [~,imx_colmag] = max(AK_rstr_glb_colmags,[],2)

            AK_raw = pagemtimes(Atns_cmp,K_glb);
            absAK = abs(AK_raw);
            absAK_evl = reshape(sum(sum(absAK,1),2),ncrv,1);
            % absAK_evl = reshape(max(sum(absAK,2),[],1),ncrv,1);
            % absAK_evl = reshape(median(sum(absAK,2),1),ncrv,1);
            % absAK_evl = reshape(median(median(absAK,2),1),ncrv,1);
            % absAK_evl = reshape(median(max(absAK,[],2),1),ncrv,1);
            % absAK_evl = reshape(max(max(absAK,[],1),[],2),ncrv,1);
            % absAK_evl = reshape(max(median(absAK,1),[],2),ncrv,1);
            % absAK_evl = reshape(median(reshape(absAK,[],ncrv),1),ncrv,1);
            [~,isort_absAK_evl] = sort(absAK_evl);
            % isort_absAK_evl = isort_absAK_evl(2:end);

            nk_use = size(K_glb,2)

            % i_sort_use = isort_absAK_evl;
            i_sort_use = isort_absAK_evl(2:end);

            [icrv_evl_mn,icrv_evl_md,icrv_evl_mx] = deal(   absAK_evl(i_sort_use(1)), ...
                                                            absAK_evl(i_sort_use(ceil(ncrv/2))), ...
                                                            absAK_evl(i_sort_use(end)))
            icrv_check = i_sort_use(1)

            icrv_evl_check = absAK_evl(icrv_check)


            % k_glb_chk = V_glb(:,end);
            % Ak_glb_chk_raw = pagemtimes(Atns_full,k_glb_chk);
            % absAk_glb_chk = abs(Ak_glb_chk_raw);
            % absAk_glb_chk_evl = reshape(sum(absAk_glb_chk,1),ncrv,1);
            % [~,isort_absAk_glb_chk_evl] = sort(absAk_glb_chk_evl);
            %
            % icrv_check2 = isort_absAk_glb_chk_evl(1);


            mspc = 'o';
            ms = 5;
            lspc = 'none';
            lw = 1;
            alpha = 0.15;

            plot(axs(1),1:ncrv,absAK_evl, ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0]);

            % absAV_plt_mat = reshape(median(abs(AV_raw),1),size(AV_raw,2),ncrv);
            absAV_plt_mat = reshape(max(abs(AV_raw),[],1),size(AV_raw,2),ncrv);
            plot(axs(2),1:size(absAV_plt_mat,1),absAV_plt_mat, ...
            'Marker','none','MarkerSize',ms,'LineStyle','-','LineWidth',1,'Color',[0 0 0 alpha]);
            plot(axs(2),1:size(absAV_plt_mat,1),absAV_plt_mat(:,icrv_), ...
            'Marker','none','MarkerSize',ms,'LineStyle','-','LineWidth',1,'Color',[nice_green 1]);
            plot(axs(2),1:size(absAV_plt_mat,1),absAV_plt_mat(:,i_sort_use(1)), ...
            'Marker','none','MarkerSize',ms,'LineStyle','-','LineWidth',1,'Color',[nice_purple 1]);

            set(axs(1),'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);
            set(axs(2),'YScale', 'linear','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);

            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,icrv_);

            spc.color = nice_green;
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_sort_use(1));

            spc.color = [0 0 1];
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_sort_use(ceil(ncrv/2)));

            spc.color = [1 0 0];
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_sort_use(end));

            spc.color = nice_purple;
            solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,icrv_check);

            solspc_plot_init_lims = solspc_plot_.get_axis_lims;

            pts_curve_check = pts_cell{icrv_check};
            vcurve_check = fspc_.ds_de_ratio_ptsmat_stabilized_fast(pts_curve_check,K_glb);
            scurve_check = [pts_curve_check(1:nvar,:) ; vcurve_check(2:(end-ndep),:)];

            spc.color = nice_pink;
            solspc_plot_ = LD_plots.plot_pts(scurve_check,meta_,solspc_plot_,spc);

            % spc.color = nice_green;
            % solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,icrv_check2);
            %
            % pts_curve_check2 = pts_cell{icrv_check2};
            % vcurve_check2 = fspc_.ds_de_ptsmat(pts_curve_check2,k_glb_chk);
            % vcurve_check2_mags = sqrt(sum(vcurve_check2.*vcurve_check2,1));
            % [vcrv_chk2_mnmag,vcrv_chk2_mdmag,vcrv_chk2_mxmag] = deal(   min(vcurve_check2_mags), ...
            %                                                             median(vcurve_check2_mags), ...
            %                                                             max(vcurve_check2_mags))
            % vcurve_check2 = vcurve_check2./vcurve_check2_mags;
            % v_scl = 1.0;
            % % v_scl = 1e8;
            % scurve_check2 = reshape([pts_curve_check2(:); pts_curve_check2(:) + v_scl*vcurve_check2(:)], ...
            %                             size(pts_curve_check2,1),2,size(pts_curve_check2,2));
            % scurve_cell_check2 = num2cell(scurve_check2,[1 2]);
            %
            % spc.color = hard_green;
            % spc.lw = 1;
            % solspc_plot_ = LD_plots.plot_pts(scurve_cell_check2,meta_,solspc_plot_,spc);
            %
            % spc.color = nice_pink;


            solspc_plot_.set_axis_lims(solspc_plot_init_lims);

            set(axs,'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);

        end
        function plt = plot_nullity_icrv_comp(S_,svd_cmp_,icrv_,solspc_plot_,plt_, cmat_)
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
            xvec = pts_mat_i(1,:);
            [npts,ncrv] = deal(size(pts_mat_i,2),S_.ncrv);
            [nrows,ncols,rho] = deal(svd_cmp_.nrow/ncrv,svd_cmp_.ncol,max(svd_cmp_.rvec));
            kappa = ncols - rho;
            nconstr_dim = nrows/npts;
            Ki = svd_cmp_.Vtns(:,(rho+1):end,icrv_);
            ATtns = reshape(svd_cmp_.matT,ncols,nrows,ncrv);

            crv_inds = 1:ncrv;
            ncrv_inds = crv_inds ~= icrv_;

            inner_K_mat = nan(ncrv,npts);
            inner_K_net = nan(ncrv,1);
            for j = 1:ncrv
                Aj = (ATtns(:,:,j))';
                % Aj = Aj./(sqrt(sum(Aj.*Aj,2)));
                abs_AjKi = abs(Aj*Ki);
                inner_K_j_sumpts = reshape(sum(reshape(abs_AjKi,nconstr_dim,npts,kappa),1),npts,kappa);

                inner_K_mat(j,:) = sum(inner_K_j_sumpts,2);
                inner_K_net(j) = sum(inner_K_mat(j,:),2);
            end
            [~,i_sort_inn] = sort(inner_K_net);

            i_mininn_net = i_sort_inn(2);

            hard_red = LD_plots.red5;
            hard_green = LD_plots.green5;
            hard_blue = LD_plots.blue5;
            nice_green = LD_plots.green4;
            nice_purple = LD_plots.purple5;
            nice_orange = LD_plots.orange1;

            blue_mat = LD_plots.get_blue_mat;
            green_mat = LD_plots.get_green_mat;
            orange_mat = LD_plots.get_orange_mat;

            cmat = cmat_;
            % cmat = blue_mat;

            spc.lw = 1.0;

            ncrv_check = 5;

            for i = 1:ncrv_check
                spc.color = cmat(i,:);
                solspc_plot_ = LD_plots.plot_solspc(S_,solspc_plot_,spc,i_sort_inn(i+1));
            end

            mspc = 'none';
            ms = 1;
            lspc = '-';
            lw = 1;
            alpha = 0.15;


            plot(axs(1),xvec,cumsum(inner_K_mat(ncrv_inds,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0 alpha]);
            for i = 1:ncrv_check
                plot(axs(1),xvec,cumsum(inner_K_mat(i_sort_inn(i+1),:),2), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',cmat(i,:));
            end

            plot(axs(2),xvec,inner_K_mat(ncrv_inds,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0 alpha]);
            for i = 1:ncrv_check
                plot(axs(2),xvec,inner_K_mat(i_sort_inn(i+1),:), ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',cmat(i,:));
            end

            % xlabel(axs(1),['$$ \sum \mathrm{err} ( \mathbf{C}^j | \mathbf{C}^i ) $$'], 'Interpreter','Latex','FontSize',14);
            % xlabel(axs(4),['$$ \mathrm{err} (s_0^j | s_0^i) $$'], 'Interpreter','Latex','FontSize',14);
            % xlabel([axs(2:3) axs(5:6)],['$$ x $$'], 'Interpreter','Latex','FontSize',14);
            xlabel(axs(1:2),['$$ x $$'], 'Interpreter','Latex','FontSize',14);

            % ylabel([axs(1) axs(4)],['$$ \sum | \mathbf{A}^j \mathbf{K}^i | $$'], 'Interpreter','Latex','FontSize',14);
            % ylabel(axs(2),['$$ \int \mathrm{err} (s^j | s^i) dx $$'], 'Interpreter','Latex','FontSize',14);
            % ylabel(axs(3),['$$ \mathrm{err} (s^j | s^i) $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(1),['$$ \int \sum | \mathbf{A}^j \mathbf{K}^i | dx $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(2),['$$ \sum | \mathbf{A}^j \mathbf{K}^i | $$'], 'Interpreter','Latex','FontSize',14);


            % set(axs(1),'YScale', 'log','XScale','log','TickLabelInterpreter','Latex','FontSize',12);
            % set([axs(1) axs(4)],'YScale', 'linear','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);
            % set([axs(2:3) axs(5:6)],'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);
            set(axs(1:2),'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);
        end
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
            alpha = 0.15;

            plot(axs(1),rel_div_net(ncrv_inds),inner_K_net(ncrv_inds), ...
            'Marker','o','MarkerSize',2,'LineStyle','none','LineWidth',lw,'Color',hard_red);
            plot(axs(1),rel_div_net(i_mindiv_ics),inner_K_net(i_mindiv_ics), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',hard_blue);
            plot(axs(1),rel_div_net(i_mindiv_net),inner_K_net(i_mindiv_net), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',hard_green);
            plot(axs(1),rel_div_net(i_mininn_net),inner_K_net(i_mininn_net), ...
            'Marker','.','MarkerSize',20,'LineStyle','none','LineWidth',lw,'Color',nice_green);

            plot(axs(2),xvec,cumsum(rel_div_mat(ncrv_inds,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0 alpha]);
            plot(axs(2),xvec,cumsum(rel_div_mat(i_mindiv_ics,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_blue);
            plot(axs(2),xvec,cumsum(rel_div_mat(i_mindiv_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_green);
            plot(axs(2),xvec,cumsum(rel_div_mat(i_mininn_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',nice_green);

            plot(axs(3),xvec,rel_div_mat(ncrv_inds,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0 alpha]);
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
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0 alpha]);
            plot(axs(5),xvec,cumsum(inner_K_mat(i_mindiv_ics,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_blue);
            plot(axs(5),xvec,cumsum(inner_K_mat(i_mindiv_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',hard_green);
            plot(axs(5),xvec,cumsum(inner_K_mat(i_mininn_net,:),2), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',nice_green);

            plot(axs(6),xvec,inner_K_mat(ncrv_inds,:), ...
            'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',[0 0 0 alpha]);
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
        function plt = plot_global_svds(pckgs_,plt_)
            plt = plt_.init_tiles_safe(1,4);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            nset = length(pckgs_);

            mspc = 'none';
            ms = 1;
            lspc = '-';
            lw = 1;
            colors = turbo(nset);

            mat_stats = @(mat_) deal(min(mat_,[],2),median(mat_,2),max(mat_,[],2), 1:size(mat_,1));

            [inds_cell,sglb_cell,inmn_cell,inmd_cell,inmx_cell] = deal(cell(nset,1));
            for i = 1:nset
                [si,Vi] = deal(pckgs_(i).svec_glb,pckgs_(i).Vmat_glb);
                sglb_cell{i} = si;
                % sglb_cell{i} = si/si(1);
                [inmn_cell{i},inmd_cell{i},inmx_cell{i},inds_cell{i}] = mat_stats(abs(pckgs_(i).matT'*Vi)');
            end


            for i = 1:nset
                leg1(i) = plot(axs(1),inds_cell{i},sglb_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end
            for i = 1:nset
                leg2(i) = plot(axs(2),inds_cell{i},inmn_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end
            for i = 1:nset
                leg3(i) = plot(axs(3),inds_cell{i},inmd_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end
            for i = 1:nset
                leg4(i) = plot(axs(4),inds_cell{i},inmx_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end

            xlabel(axs(1:4),['$$ i $$'], 'Interpreter','Latex','FontSize',14);

            ylabel(axs(1),['$$ \sigma^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(2),['min. $$ \mathbf{A}_g \mathbf{V}^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(3),['med. $$ \mathbf{A}_g \mathbf{V}^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(4),['max $$ \mathbf{A}_g \mathbf{V}^{\mathbf{A}_g}_i $$'], 'Interpreter','Latex','FontSize',14);

            % legend(axs(1), leg1,'Location', 'SouthEastOutside', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));
            legend(axs(2), leg1,'Location', 'SouthOutside', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));

            set(axs,'YScale', 'log','XScale','linear','TickLabelInterpreter','Latex','FontSize',12);

        end
        function plt = plot_curve_svds(pckgs_,plt_)
            plt = plt_.init_tiles_safe(1,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;

            nset = length(pckgs_);

            mspc = 'none';
            ms = 1;
            lspc = '-';
            lw = 1;
            colors = turbo(nset);

            mat_stats = @(mat_) deal(min(mat_,[],2),median(mat_,2),max(mat_,[],2), 1:size(mat_,1));

            [smin_cell,smed_cell,smax_cell,inds_cell] = deal(cell(nset,1));
            for i = 1:nset
                [smin_cell{i},smed_cell{i},smax_cell{i},inds_cell{i}] = mat_stats(pckgs_(i).smat);
                % [smin_cell{i},smed_cell{i},smax_cell{i},inds_cell{i}] = mat_stats(pckgs_(i).smat ./ pckgs_(i).smat(1,:) );
                mat_i = pckgs_(i).matT';
            end


            for i = 1:nset
                leg1(i) = plot(axs(1),inds_cell{i},smin_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end
            for i = 1:nset
                leg2(i) = plot(axs(2),inds_cell{i},smed_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end
            for i = 1:nset
                leg3(i) = plot(axs(3),inds_cell{i},smax_cell{i}, ...
                'Marker',mspc,'MarkerSize',ms,'LineStyle',lspc,'LineWidth',lw,'Color',colors(i,:), ...
                'DisplayName', fixlabel(pckgs_(i).dat_name));
            end

            xlabel(axs(1:3),['$$ i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(1),['min $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(2),['med. $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(3),['max $$ \sigma^{\mathbf{A}}_i $$'], 'Interpreter','Latex','FontSize',14);

            legend(axs(2), leg1,'Location', 'SouthOutside', 'Interpreter', 'Latex', 'NumColumns',min([nset,4]));

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
        function plt = plot_n1q2_solspc(pts_cell_,plt_,spc_)
            plt = plt_.init_tiles_safe(2,3);
            hold(plt.axs, 'on');
            box(plt.axs,'on');

            axs = plt.axs;

            dnames = {'x'; 'u_1'; 'u_2'; 'd_x u_1'; 'd_x u_2'};
            dim_order = [   1 2; ...
                            1 4; ...
                            2 4; ...
                            1 3; ...
                            1 5; ...
                            3 5];

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
        function plt = plot_n1q1_solspc(pts_cell_,plt_,spc_)
            plt = plt_.init_tiles_safe(1,4);
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

            view(axs(4), LD_plots.view_mat(6, :));
            for i = 1:ncrv
                plot3(axs(4), pts_cell_{i}(1,:),pts_cell_{i}(2,:), pts_cell_{i}(3,:), ...
                'Marker',spc_.mspec,'MarkerSize',spc_.ms,'LineStyle',spc_.lspec,'LineWidth',spc_.lw,'Color',spc_.color);
            end
            xlabel(axs(4),['$$' dnames{1} '$$'], 'Interpreter','Latex','FontSize',14);
            ylabel(axs(4),['$$' dnames{2} '$$'], 'Interpreter','Latex','FontSize',14);
            zlabel(axs(4),['$$' dnames{3} '$$'], 'Interpreter','Latex','FontSize',14);
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
        function plt = plot_pts(pts_,meta_,plt_,spc_)
            if (nargin==3)
                spc = LD_plots.make_default_plot_specs;
            else
                spc = spc_;
            end
            if (iscell(pts_))
                pts_cell_plot = pts_;
            else
                pts_cell_plot = mat2cell(pts_,size(pts_,1),size(pts_,2));
            end
            switch meta_.eor
                case 1
                    switch meta_.ndep
                        case 1
                            plt = LD_plots.plot_n1q1_solspc(pts_cell_plot,plt_,spc);
                        case 2
                            plt = LD_plots.plot_n1q2_solspc(pts_cell_plot,plt_,spc);
                    end
                case 2
                    switch meta_.ndep
                        case 1
                            plt = LD_plots.plot_n2q1_solspc(pts_cell_plot,plt_,spc);
                    end
            end
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
                        case 2
                            plt = LD_plots.plot_n1q2_solspc(pts_cell_plot,plt_,spc);
                    end
                case 2
                    switch S_.ndep
                        case 1
                            plt = LD_plots.plot_n2q1_solspc(pts_cell_plot,plt_,spc);
                    end
            end
        end
        function out_ = greys()
            out_ = [    LD_plots.grey1; ...
                        LD_plots.grey2; ...
                        LD_plots.grey3; ...
                        LD_plots.grey4; ...
                        LD_plots.grey5];
        end
        function out_ = purples()
            out_ = [    LD_plots.purple1; ...
                        LD_plots.purple2; ...
                        LD_plots.purple3; ...
                        LD_plots.purple4; ...
                        LD_plots.purple5];
        end
        function out_ = oranges()
            out_ = [    LD_plots.orange1; ...
                        LD_plots.orange2; ...
                        LD_plots.orange3; ...
                        LD_plots.orange4; ...
                        LD_plots.orange5];
        end
        function out_ = greens()
            out_ = [    LD_plots.green1; ...
                        LD_plots.green2; ...
                        LD_plots.green3; ...
                        LD_plots.green4; ...
                        LD_plots.green5];
        end
        function out_ = blues()
            out_ = [    LD_plots.blue1; ...
                        LD_plots.blue2; ...
                        LD_plots.blue3; ...
                        LD_plots.blue4; ...
                        LD_plots.blue5];
        end
        function out_ = reds()
            out_ = [    LD_plots.red1; ...
                        LD_plots.red2; ...
                        LD_plots.red3; ...
                        LD_plots.red4; ...
                        LD_plots.red5];
        end
    end
end

function label_out = fixlabel(label_in_)
    label_out = replace(label_in_,'_',' ');
end
