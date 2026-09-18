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
                obj = apv_plots.set_screen_posdim(obj,posdim_specs_);
            elseif (nargin == 4)
                obj = apv_plots.set_screen_posdim(obj,grid_dim_,tile_dim_,origin_tile_);
            elseif (nargin == 5)
                obj = apv_plots.set_screen_posdim(obj,grid_dim_,tile_dim_,origin_tile_,screen_);
            end
        end
        function axs_mat_out = axs_mat(obj)
            [tdim1,tdim2] = deal(obj.tile.GridSize(1),obj.tile.GridSize(2));
            axs_mat_out = (reshape(obj.axs,tdim2,tdim1))';
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
            if (iscell(lims_))
                for i = 1:naxes
                    axis(obj.axs(i),lims_{i});
                end
            else
                for i = 1:naxes
                    axis(obj.axs(i),lims_);
                end
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
        % function [plt,plt_jspc,plt_xi_spc,plt_xi_sec] = plot_LDsol_model_summary(p_,mod_,d_)
        function [plt,plt_jspc] = plot_LDsol_model_summary(p_,mod_,d_)
            nobs = size(mod_.Smat,2);
            ndep = length(mod_.fspace_0.Omap_b(:))-1;
            ndim = size(mod_.Smat,1);
            kor = ((ndim-1)/ndep) - 1;
            nvar_N1 = 1 + ndep*kor; % 1 + QN, the N1 base space dimension

            plt_jspc = apv_plots('jspc_model_summary', ...
                [ndep+1 1],...
                [ndep 1],[ndep 1], ...
                1);
            plt_jspc.show_toolbar

            function leg_out = plot_tvector(axi_,sx_,su_,vx_,vu_,clr_,LS_,mrkr_,name_)
                lw_def = 2; ms_def = 3;
                if (nargout == 0)
                    plot(axi_, ...
                        [sx_ , sx_+vx_], [su_ , su_+vu_], ...
                        'LineStyle', LS_, ...
                        'LineWidth', lw_def, ...
                        'Color', clr_, ...
                        'Marker', mrkr_, ...
                        'MarkerSize', ms_def, ...
                        'MarkerFaceColor', clr_, ...
                        'MarkerEdgeColor', clr_, ...
                        'HandleVisibility','off' ...
                    );
                else
                    if (nargin==9)
                        leg_name_ = name_;
                    else
                        leg_name_ = '';
                    end
                    leg_out = plot(axi_, ...
                        [sx_ , sx_+vx_], [su_ , su_+vu_], ...
                        'LineStyle', LS_, ...
                        'LineWidth', lw_def, ...
                        'Color', clr_, ...
                        'Marker', mrkr_, ...
                        'MarkerSize', ms_def, ...
                        'MarkerFaceColor', clr_, ...
                        'MarkerEdgeColor', clr_, ...
                        'DisplayName', leg_name_ ...
                    );
                end
            end
            function leg_out = plot2D(axi_,x_,y_,LS_,clr_,mrkr_,MS_,name_)
                if (nargout == 0)
                    plot(axi_, ...
                    x_, y_, ...
                    'LineStyle', LS_, ...
                    'LineWidth', 0.5, ...
                    'Color', clr_, ...
                    'Marker', mrkr_, ...
                    'MarkerSize', MS_, ...
                    'MarkerFaceColor', clr_, ...
                    'MarkerEdgeColor', clr_ ...
                    );
                else
                    if (nargin==8)
                        leg_name_ = name_;
                    else
                        leg_name_ = '';
                    end
                    leg_out = plot(axi_, ...
                    x_, y_, ...
                    'LineStyle', LS_, ...
                    'LineWidth', 0.5, ...
                    'Color', clr_, ...
                    'Marker', mrkr_, ...
                    'MarkerSize', MS_, ...
                    'MarkerFaceColor', clr_, ...
                    'MarkerEdgeColor', clr_, ...
                    'DisplayName', leg_name_ ...
                    );
                end
            end
            function plot_labelled_svds(axi_,svds_,lbls_,leg_loc_,nrm_flg_)
                if (nargin==4)
                    nrm_flg = false;
                else
                    nrm_flg = nrm_flg_;
                end
                if (nrm_flg)
                    s_plt = @(s_) s_/(s_(1));
                else
                    s_plt = @(s_) s_;
                end
                nsvd = length(svds_(:));
                cmat = hsv(nsvd);
                if (nsvd==1)
                    leg_i(1) = plot(axi_, ...
                    1:length(svds_.s), s_plt(svds_.s), ...
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
                        1:length(svds_{i}.s), s_plt(svds_{i}.s), ...
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
                    % 'XScale', 'log', ...
                xlabel(axi_, '$$ i $$', 'Interpreter','Latex','FontSize',16);
                ylabel(axi_, '$ \sigma_i $', 'Interpreter','Latex','FontSize',16);
                legend(axi_, leg_i(1:nsvd),'Location', leg_loc_, 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',14);
            end
            rlbl = @(s_) [ ', \rho=' num2str(s_.r) '/' num2str(s_.dim) ];
            nlbl = @(s_) [ ', \| \cdot \|_* =' num2str(sum(s_.s/s_.s(1)),'%.1f') '/' num2str(s_.dim) ];

            [plt_jspc,pspc,leg_S_1] = apv_plots.plot_Sobs(plt_jspc,mod_.Sobs,d_);
            lims0 = plt_jspc.get_axis_lims();

            pspc_i = pspc;
            pspc_i.Color = [1 1 1];
            [plt_jspc,pspc_i,leg_S_2] = apv_plots.plot_Sobs(plt_jspc, mod_.Sobs{mod_.icrv_sO},pspc_i);
            pspc_ii = pspc_i;
            plt = plt_jspc;
            axs = plt.axs;
            axs_mat = plt.axs_mat;

            ileg = 1;
            for k = 1:(kor+1)
                if (k==1)
                    dkx_str_k = '';
                else
                    dkx_str_k = ['d^{' num2str(k-1) '}'];
                end
                for i = 1:ndep
                    dkxui_str_ki = [ dkx_str_k 'u_{' num2str(i) '}'];
                    leg_S_1(ileg).DisplayName = ['$ (x , ' dkxui_str_ki ')|_S $'];
                    leg_S_2(ileg).DisplayName = ['$ (x , ' dkxui_str_ki ')|_{C_{\textrm{origin}}} $'];

                    legend(plt_jspc.axs(ileg), [ ...
                        leg_S_1(ileg), ...
                        leg_S_2(ileg) ...
                        ],'Interpreter', 'Latex', ...
                        'Location', 'EastOutside', ...
                        'NumColumns',1,'FontSize',12);

                    ileg = ileg + 1;
                end
            end
            dkx_str = @(k_) ['d_x^{' num2str(k_) '}'];
            ui_str = @(i_) ['u_{' num2str(i_) '}'];
            dkxui_str = @(k_,i_) [ dkx_str(k_) ui_str(i_) ];
            xii_str = @(i_) ['\xi_{' num2str(i_) '}'];
            upy_str = @(y_) ['\upsilon^{' num2str(y_) '}'];

            tau_uN_RN1_tns = reshape(mod_.tau_uN_RN1_net(:,1,:),ndep,kor,[]);
            sO = mod_.s_O;
            [xO,uO] = deal( mod_.s_O(1),reshape(mod_.sNp1_O(2:end),ndep,kor+2) );
            t_sO = mod_.t_O;

            % M_sO_basis = mod_.Mcom_sO_basis;
            % M_sO_basis = mod_.Mnet_sO_basis;

            % A_G_sO = M_sO_basis.U * diag(M_sO_basis.s) * M_sO_basis.V';
            % A_G_sO = mod_.flow_pckg. .* ()';
            % A_G_sO = mod_.flow_pckg.Tspc_N_sO_svds(1) .* ()';
            % A_G_sO = A_G_sO * ( norm(t_sO(1:ndim)) / sqrt(max(sum(A_G_sO.*A_G_sO,1))) ); % rescale wrt tvf

            % A_G_sO = (mod_.Gnet_sO_coords.s') .* mod_.Gnet_sO_coords.U;
            % A_G_sO = mod_.Gcom_sO_coords.sO_nTVF.Tspc_image;
            % A_G_sO = mod_.Gcom_sO_coords.sO_nTVF.nTVF_Tspc_image;
            % A_G_sO = mod_.Gnet_sO_coords.sO_nTVF.Tspc_image;
            % A_G_sO = mod_.Gnet_sO_coords.sO_nTVF.nTVF_Tspc_image;

            % A_G_sO = mod_.Ncom_sNO_basis.sO_nTVF.Tspc_image;
            % A_G_sO = mod_.Ncom_sNO_basis.sO_nTVF.nTVF_Tspc_image;
            % A_G_sO = mod_.Nnet_sNO_basis.sO_nTVF.Tspc_image;
            % A_G_sO = mod_.Nnet_sNO_basis.sO_nTVF.nTVF_Tspc_image;

            % A_G_sO = mod_.flow_pckg.VN_spc_sO;
            % A_G_sO = mod_.flow_pckg.VN_spc_sO(:,2:end);
            A_G_sO = mod_.flow_pckg.Gn_bse.Tspc_Nv_sO_mat0;
            % A_G_sO = mod_.flow_pckg.Tspc_Nv_sO_tns(:,:,1);
            % A_G_sO = A_G_sO./sqrt(sum(A_G_sO.^2,1));

            % A_G_sO = nan(ndim,length(mod_.flow_pckg.Tspc_N_sO_svds(:)));
            % for i = 1:size(A_G_sO,2)
            %     A_G_sO(:,i) = sum(mod_.flow_pckg.Tspc_N_sO_svds(i).U .*  mod_.flow_pckg.Tspc_N_sO_svds(i).s' , 2);
            %     A_G_sO(:,i) = A_G_sO(:,i)/norm(A_G_sO(:,i));
            % end

            A_G_sO = A_G_sO * ( norm(t_sO(1:ndim)) / sqrt(max(sum(A_G_sO.*A_G_sO,1))) ); % rescale wrt tvf
            A_G_sO_u = reshape(A_G_sO(2:end,:),ndep,kor+1,[]);

            % Aspc.cmat = apv_plots.orange1 .* ones(size(A_G_sO,2),3);
            % Aspc.cmat = nebula(size(A_G_sO,2));
            Aspc.cmat = autumn(size(A_G_sO,2));
            % Aspc.cmat = cool(size(A_G_sO,2));
            % Aspc.cmat = spring(size(A_G_sO,2));
            Aspc.lspc = ':';
            Aspc.mspc = 'd';

            % Vspc_O = M_sO_basis.Vspc_sO;
            % Vspc_O = mod_.flow_pckg.VN_spc_sO;
            Vspc_O = mod_.flow_pckg.VN_spc_sO ./ sqrt(sum(mod_.flow_pckg.VN_spc_sO.^2,1));

            % Vspc_svd_i = mod_.flow_pckg.Tspc_N_sO_svds(1);
            % Vspc_O = M_sO_basis.Vspc_sO;
            % Vspc_svd_i = mod_.flow_pckg.Gn_bse.TTspc_svd_0;
            % Vspc_O = Vspc_svd_i.U(:,1:nvar_N1) .* ( Vspc_svd_i.s(1:nvar_N1) )';

            Vspc_O = Vspc_O*( norm(t_sO(1:ndim)) / sqrt(max(sum(Vspc_O.*Vspc_O,1))) ); % rescale wrt tvf
            Vspc_O_x = Vspc_O(1,:);
            Vspc_O_u = reshape( Vspc_O(2:end,:), ndep,kor+1,[] );
            % J_xi_sO = [ M_sO_basis.JXi_sO_TspcO ; zeros( ndep,nvar_N1 ) ];
            % J_xi_sO = J_xi_sO*( norm(t_sO(1:ndim)) / sqrt(max(sum(J_xi_sO.*J_xi_sO,1))) ); % rescale wrt tvf
            % J_xi_O_x = J_xi_sO(1,:);
            % J_xi_O_u = reshape( J_xi_sO(2:end,:), ndep,kor+1,[] );
            % cmat_i = hsv(nvar_N1);
            cmat_i = autumn(nvar_N1);
            for i = 1:ndep
                tau_uiN = reshape(tau_uN_RN1_tns(i,:,:),kor,[]);

                ui_str_i = [ 'u_{' num2str(i) '}'];
                for itspc = 1:size(A_G_sO,2)
                    plot_tvector( ...
                        axs_mat(i,1),xO,uO(i,1),A_G_sO(1,itspc),A_G_sO_u(i,1,itspc),Aspc.cmat(itspc,:),Aspc.lspc,Aspc.mspc);
                end
                leg_ii(i,1,1) = plot2D( ...
                    axs_mat(i,1),xO,uO(i,1),'none',[1 1 1],'o',pspc.MarkerSize, ...
                    ['$ ( x, ' ui_str_i ')|_{s_O} $'] ...
                );
                % ['$ ( x, \tau_{' dkxui_str(k,i) '} ) |_S := ( x, \Lambda_{' dkxui_str(k,i) '} \vartheta )|_S $
                leg_ii(i,2,1) = plot_tvector( ...
                    axs_mat(i,1),xO,uO(i,1),1,uO(i,2),[1 1 1],'-','none', ...
                    ['$ (1, \tau_{' ui_str_i '} ) |_{s_O} ) $'] ...
                );
                ileg = 2;
                for iq = 1:nvar_N1
                    % ['$ ( ' upy_str(iq) '_x,' upy_str(iq) '_{' dkxui_str(k-1,i) '} ) |_{s_O} := d_{' xii_str(iq) '} (x,' dkxui_str(k-1,i) ') |_{s_O} $'] ...
                    leg_ii(i,ileg+1,1) = plot_tvector( ...
                        axs_mat(i,1),xO,uO(i,1),Vspc_O_x(iq),Vspc_O_u(i,1,iq),cmat_i(iq,:),'-','none', ...
                        ['$ ( ' upy_str(iq) '_x,' upy_str(iq) '_{' ui_str_i '} ) |_{s_O} $'] ...
                    );
                    ileg = ileg + 1;
                    % leg_ii(i,ileg+2,1) = plot_tvector( ...
                    %     axs_mat(i,1),xO,uO(i,1),J_xi_O_x(iq),J_xi_O_u(i,1,iq),cmat_i(iq,:),':','d', ...
                    %     ['$ ( \partial_x' xii_str(iq) ', \partial_{' ui_str_i '}' xii_str(iq) ') |_{s_O} $'] ...
                    % );
                    % ileg = ileg + 2;
                end

                % ['$ ( x, \tau_{' dkxui_str(k,i) '} ) |_S := ( x, \Lambda_{' dkxui_str(k,i) '} \vartheta )|_S $
                for k = 2:(kor+1)
                    % axs_mat(i,k),mod_.Smat(1,:),tau_uiN(k-1,:),'none',[0 1 0],'o',ceil(0.25*pspc.MarkerSize), ...
                    leg_ii(i,1,k) = plot2D( ...
                        axs_mat(i,k),mod_.Smat(1,:),tau_uiN(k-1,:),'none',[0 1 0],'o',ceil(0.25*pspc.MarkerSize), ...
                        ['$ ( x, \Lambda_{' dkxui_str(k,i) '} \vartheta )|_S $ '] ...
                    );
                    for itspc = 1:size(A_G_sO,2)
                        plot_tvector( ...
                            axs_mat(i,k),xO,uO(i,k),A_G_sO(1,itspc),A_G_sO_u(i,k,itspc),Aspc.cmat(itspc,:),Aspc.lspc,Aspc.mspc);
                    end
                    leg_ii(i,2,1) = plot2D( ...
                        axs_mat(i,k),xO,uO(i,k),'none',[1 1 1],'o',pspc.MarkerSize, ...
                        ['$ ( x, ' dkxui_str(k-1,i) ')|_{s_O} $'] ...
                    );
                    % ['$ ( \tau_x, \tau_{' dkxui_str(k-1,i) '} ) |_{s_O} := (1, ' dkxui_str(k,i) ' |_{s_O} ) $'] ...
                    leg_ii(i,3,k) = plot_tvector( ...
                        axs_mat(i,k),xO,uO(i,k),1,uO(i,k+1),[1 1 1],'-','none', ...
                        ['$ (1, \tau_{' dkxui_str(k-1,i) '} ) |_{s_O} ) $'] ...
                    );
                    ileg = 3;
                    for iq = 1:nvar_N1
                        % ['$ ( ' upy_str(iq) '_x,' upy_str(iq) '_{' dkxui_str(k-1,i) '} ) |_{s_O} := d_{' xii_str(iq) '} (x,' dkxui_str(k-1,i) ') |_{s_O} $'] ...
                        leg_ii(i,ileg+1,k) = plot_tvector( ...
                            axs_mat(i,k),xO,uO(i,k),Vspc_O_x(iq),Vspc_O_u(i,k,iq),cmat_i(iq,:),'-','none', ....
                            ['$ ( ' upy_str(iq) '_x,' upy_str(iq) '_{' dkxui_str(k-1,i) '} ) |_{s_O} $'] ...
                        );
                        ileg = ileg + 1;
                        % leg_ii(i,ileg+2,k) = plot_tvector( ...
                        %     axs_mat(i,k),xO,uO(i,k),J_xi_O_x(iq),J_xi_O_u(i,k,iq),cmat_i(iq,:),':','d', ...
                        %     ['$ ( \partial_x' xii_str(iq) ', \partial_{' dkxui_str(k-1,i) '}' xii_str(iq) ') |_{s_O} $'] ...
                        % );
                        % ileg = ileg + 2;
                    end
                end
            end
            if ( (ndep==kor)&&(kor==1) )
                sO_11 = reshape(sO(1:3),[],1);
                plt3D_i = @(s_,mrkr_,c_,LS_) plot3(axs_mat(1,end), ...
                    s_(1,:), s_(2,:), s_(3,:), ...
                    'LineStyle', LS_, ...
                    'LineWidth', 1, ...
                    'Color', c_, ...
                    'Marker', mrkr_, ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', c_, ...
                    'MarkerEdgeColor', [0 0 0], ...
                    'HandleVisibility','off' ...
                );
                for itspc = 1:size(A_G_sO,2)
                    plt3D_i([ sO_11, sO_11+[ A_G_sO(1,itspc) ; A_G_sO_u(1,1:2,itspc)' ] ],'none',Aspc.cmat(itspc,:),Aspc.lspc);
                end
                plt3D_i(sO_11,'o',[1 1 1],'none')
                plt3D_i([ sO_11, sO_11+[ t_sO(1) ; t_sO(2:3) ] ],'d',[1 1 1],'-')
                plt3D_i([ sO_11, sO_11+[ Vspc_O_x(1) ; Vspc_O_u(1,1:2,1)' ] ],'d',cmat_i(1,:),'-')
                plt3D_i([ sO_11, sO_11+[ Vspc_O_x(2) ; Vspc_O_u(1,1:2,2)' ] ],'d',cmat_i(2,:),'-')
                % plt3D_i([ sO_11, sO_11+[ J_xi_O_x(1) ; J_xi_O_u(1,1:2,1)' ] ],'d',cmat_i(1,:),':')
                % plt3D_i([ sO_11, sO_11+[ J_xi_O_x(2) ; J_xi_O_u(1,1:2,2)' ] ],'d',cmat_i(2,:),':')
            end
            % keyboard

            %% visualize model diagnostics
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
            % plt.show_menubar
            plt.show_toolbar

            axi = axs_mat(1,1);
            title(axi, ['SVDs over parameter space '], ...
                'Interpreter','Latex','FontSize',12 );
            svds_i = { ...
                mod_.Jl_N1_svd , '$J_{\lambda} |_{ \{ s_j \} } '; ...
                mod_.H_N1_svd , '$H^{(0)}_{\textrm{tvf}}'; ...
                mod_.DprN_svd, '$D_{\textrm{pr}}'; ...
                mod_.Rsvd_N1_net, '$R^1_{\mathrm{net}}'; ...
                mod_.Gsvd_N1_net,  '$G_{\mathrm{net}}'; ...
                mod_.Bsvd_t0_N1,  '$B^{(0)}'; ...
                mod_.Bsvd_tdxu_N1,  '$B_{d_x u}'; ...
                mod_.Gsvd_N1_com,  '$G_{\mathrm{com}}'; ...
                mod_.Bsvd_tN_N1,  '$B^{(1)}'; ...
                mod_.Nsvd_N1_net,  '$N_{\mathrm{net}}'; ...
                mod_.Nsvd_N1_com,  '$N_{\mathrm{com}}'; ...
            };
            nsvd_i = size(svds_i,1);
            Bsvds_i = mod_.flow_pckg.B_v_svds;
            for isvd = 1:length(Bsvds_i(:))
                svds_i{isvd+nsvd_i,1} = Bsvds_i(isvd);
                svds_i{isvd+nsvd_i,2} = ['$B_{' num2str(isvd) '}'];
            end
            nsvd_i = size(svds_i,1);
            Nsvds_i = mod_.flow_pckg.N_v_svds;
            for isvd = 1:length(Nsvds_i(:))
                svds_i{isvd+nsvd_i,1} = Nsvds_i(isvd);
                svds_i{isvd+nsvd_i,2} = ['$N_{' num2str(isvd) '}'];
            end
            nsvd_i = size(svds_i,1);
            Gsvds_i = mod_.flow_pckg.Gc_v_svds;
            for isvd = 1:length(Gsvds_i(:))
                svds_i{isvd+nsvd_i,1} = Gsvds_i(isvd);
                svds_i{isvd+nsvd_i,2} = ['$G_c^{' num2str(isvd) '}'];
            end
            nsvd_i = size(svds_i,1);

            svds = cell([ size(svds_i,1),1 ]);
            labels = cell([length(svds),1]);
            for isvd = 1:size(svds_i,1)
                labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) nlbl(svds_i{isvd,1}) '$'];
                svds{isvd} = svds_i{isvd,1};
            end
            plot_labelled_svds(axi,svds,labels,'EastOutside',true);

            axi = axs_mat(1,2);
            title(axi, ['SVDs over base/jet space (at the origin) '], ...
                'Interpreter','Latex','FontSize',12 );
            svds_i = { ...
                mod_.Jl_N1_svd.gEta_tvf_YJl_sO_svd, '$\nabla^{(0)} ( \eta_i )_{i=1}^{\rho} '; ...
                mod_.Mcom_sO_basis , '$\Lambda^{(1)}_{O} W_{G_c}'; ...
                mod_.Mnet_sO_basis , '$\Lambda^{(1)}_{O} W_{G_n}'; ...
                mod_.Ncom_sO_basis , '$\Lambda^{(1)}_{O} W_{G_c} D_{N_c}'; ...
                mod_.Nnet_sO_basis , '$\Lambda^{(1)}_{O} W_{G_n} D_{N_n}'; ...
            };
            nsvd_i = size(svds_i,1);

            Tsvds_i = mod_.flow_pckg.Tspc_N_sO_svds;
            for isvd = 1:length(Tsvds_i(:))
                svds_i{isvd+nsvd_i,1} = Tsvds_i(isvd);
                svds_i{isvd+nsvd_i,2} = ['$T_{' num2str(isvd) '}'];
            end
            nsvd_i = size(svds_i,1);

            svds = cell([ size(svds_i,1),1 ]);
            labels = cell([length(svds),1]);
            for isvd = 1:size(svds_i,1)
                labels{isvd} = [ svds_i{isvd,2} rlbl(svds_i{isvd,1}) nlbl(svds_i{isvd,1}) '$'];
                svds{isvd} = svds_i{isvd,1};
            end
            plot_labelled_svds(axi,svds,labels,'EastOutside',true);


            axi = axs_mat(2,1);
            title(axi, ['trivial vector field error'], ...
                'Interpreter','Latex','FontSize',12 );
            cmat_i = hsv(kor*ndep);
            ileg = 1;
            for k = 1:kor
                for i = 1:ndep
                    leg_i(ileg) = plot(axi, ...
                        mod_.Smat(1,:) , reshape(mod_.err_tau_u_tvf(i,k,:),1,[]), ...
                        'LineStyle', 'none', ...
                        'LineWidth', 0.5, ...
                        'Color', cmat_i(ileg,:), ...
                        'Marker', 'o', ...
                        'MarkerSize', pspc.MarkerSize, ...
                        'MarkerFaceColor', cmat_i(ileg,:), ...
                        'MarkerEdgeColor', cmat_i(ileg,:), ...
                        'DisplayName', ['$ d^' num2str(k) '_x u_' num2str(i) '$'] ...
                    );
                    ileg = ileg+1;
                end
            end
            set(axi, ...
                'XScale', 'linear' , ...
                'YScale', 'linear' );
            xlabel(axi, '$ x $', 'Interpreter','Latex','FontSize',16);
            ylabel(axi, '$ d^k_x u_i $ error (tvf model - observed)', 'Interpreter','Latex','FontSize',16);
            legend(axi, leg_i(1:(kor*ndep)),'Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',14);

            axi = axs_mat(2,2);
            title(axi, ['trivial vector field canonical dependent coordinates, satisfying $ \mathbf{\tau} \cdot \nabla \eta_{\textrm{tvf}} = d_{x} (\eta_{\textrm{tvf}}) = 0 $'], ...
                'Interpreter','Latex','FontSize',12 );
            nleg_i = size(mod_.Jl_N1_svd.Eta_tvf_S,2);
            cmat_i = cool(nleg_i);
            for i = 1:nleg_i
                leg_i(i) = plot(axi, ...
                    mod_.Smat(1,:), ...
                    reshape( ...
                    abs(mod_.Jl_N1_svd.Eta_tvf_S(:,i) - mod_.Jl_N1_svd.Eta_tvf_sO(i)), ...
                    1, []), ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', cmat_i(i,:), ...
                    'Marker', 'o', ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', cmat_i(i,:), ...
                    'MarkerEdgeColor', cmat_i(i,:), ...
                    'DisplayName', ['$ \eta_{\textrm{tvf}}^{' num2str(i) '}$'] ...
                );
            end
            legend(axi, leg_i(1:nleg_i),'Location', 'EastOutside', 'Interpreter', 'Latex', 'NumColumns',1,'FontSize',14);
            set(axi, ...
                'XScale', 'linear' , ...
                'YScale', 'log' );
            xlabel(axi, '$ x $', 'Interpreter','Latex','FontSize',16);
            ylabel(axi, '$ \left\| \eta^i_{\textrm{tvf}} |_S - \eta^i_{\textrm{tvf}} |_{s_O} \right\| $', 'Interpreter','Latex','FontSize',16);

            %% visualize Sobs over canonical coordinates
            % if (ndim==3)
            %     plt_Xispc_width = 2;
            % else
            %     plt_Xispc_width = 1;
            % end
            % plt_xi_spc = apv_plots('Xi_spc', ...
            %     [ndim 2],...
            %     [nvar_N1 plt_Xispc_width],[ndim 1], ...
            %     1);
            % plt_xi_spc.show_toolbar
            % d_xi = pspc;
            % % d_xi.ndep = nvar_N1-1;
            % d_xi.ndep = nvar_N1;
            % d_xi.eor = 0;
            % d_xi.LineStyle = 'none';
            % Sobs_Xi = mod_.dXi_S_sO_cell;
            % for i = 1:length(mod_.Sobs)
            %     Sobs_Xi{i} = [mod_.Sobs{i}(1,:) ; Sobs_Xi{i}];
            % end
            % [plt_xi_spc,pspc_xi,leg_xi_1] = apv_plots.plot_Sobs(plt_xi_spc,Sobs_Xi,d_xi);
            % if (mod_.isrtmags_dXi_S_sO_crv(1) ~= mod_.icrv_sO)
            %     pspc_xi.Color = [1 1 1];
            %     [plt_xi_spc,pspc_xi] = apv_plots.plot_Sobs(plt_xi_spc,Sobs_Xi{mod_.icrv_sO},pspc_xi);
            % end
            % pspc_xi.Color = apv_plots.blue1;
            % [plt_xi_spc,pspc_xi,leg_xi_2] = apv_plots.plot_Sobs(plt_xi_spc,Sobs_Xi{mod_.isrtmags_dXi_S_sO_crv(1)},pspc_xi);
            % pspc_xi.Color = apv_plots.red1;
            % [plt_xi_spc,pspc_xi,leg_xi_3] = apv_plots.plot_Sobs(plt_xi_spc,Sobs_Xi{mod_.isrtmags_dXi_S_sO_crv(2)},pspc_xi);
            % pspc_xi.Color = [1 1 1];
            % [plt_xi_spc,pspc_xi,leg_xi_4] = apv_plots.plot_Sobs(plt_xi_spc,[ sO(1) ; zeros(nvar_N1,1) ],pspc_xi);
            % pspc_ii = pspc_xi;
            % pspc_ii.Color = apv_plots.grey3;
            % pspc_ii.LineStyle = 'none';
            % Splt_i = [ mod_.Smat(1,:) ; mod_.dXi_S_sO ];
            % [plt_xi_spc,~,leg_xi_5] = apv_plots.plot_Sobs(plt_xi_spc,Splt_i(:,mod_.isrtmags_dXi_S_sO(1:(2*nvar_N1))),pspc_ii);
            % for i = 1:nvar_N1
            %     xii_str_i = ['\xi_{' num2str(i) '}'];
            %     leg_xi_1(i).DisplayName = ['$ (x ,' xii_str_i '|_S-' xii_str_i '|_{s_O})$'];
            %     leg_xi_2(i).DisplayName = ['$ (x ,' xii_str_i '|_{C_{\textrm{origin}}}-' xii_str_i '|_{s_O})$'];
            %     leg_xi_3(i).DisplayName = ['$ (x ,' xii_str_i '|_{C_{\textrm{neighbor}}}-' xii_str_i '|_{s_O})$'];
            %     leg_xi_4(i).DisplayName = ['$ (x ,' xii_str_i '|_{s_O}-' xii_str_i '|_{s_O})$'];
            %     leg_xi_5(i).DisplayName = ['$ (x ,' xii_str_i '|_{\Xi_{\textrm{neighbor}}}-' xii_str_i '|_{s_O})$'];
            %     legend(plt_xi_spc.axs(i), [ leg_xi_1(i) leg_xi_2(i) leg_xi_3(i) leg_xi_4(i) leg_xi_5(i) ],'Interpreter', 'Latex', ...
            %     'Location', 'EastOutside', ...
            %     'NumColumns',1,'FontSize',14);
            %     ylabel(plt_xi_spc.axs(i), ['$ ' xii_str(i) '-' xii_str(i) '|_{s_O} $'], 'Interpreter','Latex','FontSize',16);
            % end

            % plt_xi_sec = apv_plots('Xi_spc_sections', ...
            %     [1 1],...
            %     [1 1],[1 1], ...
            %     1);
            % [ndim 1],...
            % [nvar_N1 1],[ndim 1], ...
            % fntsze_Xi_sec = 12;
            % plt = plt_xi_sec;
            % plt = plt.init_tiles_safe(nvar_N1,nvar_N1);
            % hold(plt.axs, 'on');
            % box(plt.axs,'on');
            % axs = plt.axs;
            % axs_mat = plt.axs_mat;
            % set(axs, ...
            %     'YScale', 'linear',  ...
            %     'XScale', 'linear',  ...
            %     'TickLabelInterpreter','Latex', ...
            %     'FontSize',16 );
            % plt.show_toolbar
            % % Hsvd_tvf_YHnet_S = M_sO_basis.Hsvd_tvf_YHnet_S;
            % Smat_xi = M_sO_basis.Xi_S - M_sO_basis.Xi_sO;
            % for i = 1:nvar_N1
            %     xii_str_i = ['\xi_{' num2str(i) '}'];
            %     for j = 1:(i-1)
            %         leg_i(1) = plot2D(axs_mat(i,j), ...
            %             Smat_xi(j,:),Smat_xi(i,:), ...
            %             'none', cmat_i(i,:), 'o', pspc.MarkerSize, ...
            %             ['$ ( ' xii_str(j) '|_S -' xii_str(j) '|_{s_O} ,' xii_str(i) '|_S - ' xii_str(i) ' |_{s_O} ) $'] ...
            %         );
            %         leg_i(2) = fplot(axs_mat(i,j), ...
            %             @(x_) x_ , [min(Smat_xi(j,:)),max(Smat_xi(j,:))],  ...
            %             'LineStyle', ':', ...
            %             'Color', [1 1 1], ...
            %             'DisplayName', ['$ ' xii_str(j) ' - ' xii_str(j) '|_{s_O} $'] ...
            %         );
            %         leg_i(3) = plot2D( axs_mat(i,j), ...
            %             0,0, ...
            %             'none',[1 1 1], 'o', pspc.MarkerSize, ...
            %             ['$ (0, 0) $'] ...
            %         );
            %         legend(axs_mat(i,j), leg_i(1:2),'Interpreter', 'Latex', ...
            %         'Location', 'NorthOutside', ...
            %         'NumColumns',1,'FontSize',fntsze_Xi_sec);
            %         xlabel(axs_mat(i,j), ['$ ' xii_str(j) ' - ' xii_str(j) '|_{s_O} $'], 'Interpreter','Latex','FontSize',fntsze_Xi_sec);
            %         axis(axs_mat(i,j), [min(Smat_xi(j,:)),max(Smat_xi(j,:)),min(Smat_xi(i,:)),max(Smat_xi(i,:))]);
            %     end
            %     leg_i(1) = plot2D(axs_mat(i,i), ...
            %         mod_.Smat(1,:) - sO(1),Smat_xi(i,:), ...
            %         'none', cmat_i(i,:), 'o', pspc.MarkerSize, ...
            %         ['$ ( x |_S - x |_{s_O} ,' xii_str_i '|_S -' xii_str_i ' |_{s_O} ) $'] ...
            %     );
            %     leg_i(2) = fplot(axs_mat(i,i), ...
            %         @(x_) x_ , [min(mod_.Smat(1,:) - sO(1)),max(mod_.Smat(1,:) - sO(1))],  ...
            %         'LineStyle', ':', ...
            %         'Color', [1 1 1], ...
            %         'DisplayName', ['$ x - x|_{s_O} $'] ...
            %     );
            %     leg_i(3) = plot2D( axs_mat(i,i), ...
            %         0,0, ...
            %         'none',[1 1 1], 'o', pspc.MarkerSize, ...
            %         ['$ (0, 0) $'] ...
            %     );
            %     legend(axs_mat(i,i), leg_i(1:2),'Interpreter', 'Latex', ...
            %     'Location', 'NorthOutside', ...
            %     'NumColumns',1,'FontSize',fntsze_Xi_sec);
            %     xlabel(axs_mat(i,i), ['$ x - x|_{s_O} $'],'Interpreter','Latex','FontSize',fntsze_Xi_sec);
            %     axis(axs_mat(i,i), [min(mod_.Smat(1,:)),max(mod_.Smat(1,:)),min(Smat_xi(i,:)),max(Smat_xi(i,:))]);
            %     for j = (i+1):nvar_N1
            %         leg_i(1) = plot2D(axs_mat(i,j), ...
            %             Smat_xi(j,:),Smat_xi(i,:), ...
            %             'none', cmat_i(i,:), 'o', pspc.MarkerSize, ...
            %             ['$ ( ' xii_str(j) '|_S -' xii_str(j) '|_{s_O} ,' xii_str(i) '|_S - ' xii_str(i) ' |_{s_O} ) $'] ...
            %         );
            %         leg_i(2) = fplot(axs_mat(i,j), ...
            %             @(x_) x_ , [min(Smat_xi(j,:)),max(Smat_xi(j,:))],  ...
            %             'LineStyle', ':', ...
            %             'Color', [1 1 1], ...
            %             'DisplayName', ['$ ' xii_str(j) ' - ' xii_str(j) '|_{s_O} $'] ...
            %         );
            %         leg_i(3) = plot2D( axs_mat(i,j), ...
            %             0,0, ...
            %             'none',[1 1 1], 'o', pspc.MarkerSize, ...
            %             ['$ (0, 0) $'] ...
            %         );
            %         legend(axs_mat(i,j), leg_i(1:2),'Interpreter', 'Latex', ...
            %         'Location', 'NorthOutside', ...
            %         'NumColumns',1,'FontSize',fntsze_Xi_sec);
            %         xlabel(axs_mat(i,j), ['$ ' xii_str(j) ' - ' xii_str(j) '|_{s_O} $'], 'Interpreter','Latex','FontSize',fntsze_Xi_sec);
            %         axis(axs_mat(i,j), [min(Smat_xi(j,:)),max(Smat_xi(j,:)),min(Smat_xi(i,:)),max(Smat_xi(i,:))]);
            %     end
            %     ylabel(axs_mat(i,:), ['$ ' xii_str(i) '-' xii_str(i) '|_{s_O} $'], 'Interpreter','Latex','FontSize',fntsze_Xi_sec);
            % end

            % plt_xi3 = apv_plots('Xi3', ...
            %     [3 3],...
            %     [1 1],[1 3], ...
            %     1);
            % plt = plt_xi3;
            % plt = plt.init_tiles_safe(1,1);
            % hold(plt.axs, 'on');
            % box(plt.axs,'on');
            % axs = plt.axs;
            % axs_mat = plt.axs_mat;
            % set(axs, ...
            %     'YScale', 'linear',  ...
            %     'XScale', 'linear',  ...
            %     'TickLabelInterpreter','Latex', ...
            %     'FontSize',16 );
            % plt.show_toolbar
            % plt3D_i = @(s_,mrkr_,c_,LS_) plot3(axs(1), ...
            %     s_(1,:), s_(2,:), s_(3,:), ...
            %     'LineStyle', LS_, ...
            %     'LineWidth', 1, ...
            %     'Color', c_, ...
            %     'Marker', mrkr_, ...
            %     'MarkerSize', pspc.MarkerSize, ...
            %     'MarkerFaceColor', c_, ...
            %     'MarkerEdgeColor', [0 0 0] ...
            % );
            % if (nvar_N1>=3)
            %     for i = 1:size(mod_.ipts_crv,2)
            %         plt3D_i(Smat_xi(:,mod_.ipts_crv(1,i):mod_.ipts_crv(2,i)),'o',apv_plots.green4,'-');
            %     end
            %     view(axs(1), apv_plots.view_mat(6, :));
            %     xlabel(axs(1), '$ \xi_1 $', 'Interpreter','Latex','FontSize',16);
            %     ylabel(axs(1), '$ \xi_2 $', 'Interpreter','Latex','FontSize',16);
            %     zlabel(axs(1), '$ \xi_3 $', 'Interpreter','Latex','FontSize',16);
            % else
            %     for i = 1:size(mod_.ipts_crv,2)
            %         plt3D_i( [ mod_.Smat(1,mod_.ipts_crv(1,i):mod_.ipts_crv(2,i)) ; ...
            %             Smat_xi(:,mod_.ipts_crv(1,i):mod_.ipts_crv(2,i)) ...
            %             ] , 'o',apv_plots.green4,'-');
            %     end
            %     view(axs(1), apv_plots.view_mat(6, :));
            %     xlabel(axs(1), '$ \xi_1 $', 'Interpreter','Latex','FontSize',16);
            %     ylabel(axs(1), '$ \xi_2 $', 'Interpreter','Latex','FontSize',16);
            %     zlabel(axs(1), '$ \xi_3 $', 'Interpreter','Latex','FontSize',16);
            % end

            plt_jspc.set_axis_lims(lims0);

        end
        function [plt,pspc,leg_out] = plot_Sobs(p_,S_,d_)
            plt = p_;

            pspc = apv_plots.verify_plotspecs(d_);

            eor = d_.eor;
            ndep = d_.ndep;
            ndim = 1+ndep*(eor+1);

            if (isempty(plt.axs))
                if (ndim==3)
                    plt = plt.init_tiles_safe(1,3);
                else
                    plt = plt.init_tiles_safe(ndep,eor+1);
                end
            end
            hold(plt.axs, 'on');
            box(plt.axs,'on');
            axs = plt.axs;
            axs_mat = plt.axs_mat;
            if (eor==0) % canonical coordinate plot
                xlabel(axs(1:ndep), '$$ x $$', 'Interpreter','Latex','FontSize',16);
                yname = @(i_,k_) ['$$ \xi_' num2str(i_) '$$'];
                xname_3D = '$$ \xi_1 $$';
                yname_3D = '$$ \xi_2 $$';
                zname_3D = '$$ \xi_3 $$';
            else
                axs_set = axs_mat(:,1:(eor+1));
                set(axs_set, ...
                    'YScale', 'linear',  ...
                    'XScale', 'linear',  ...
                    'TickLabelInterpreter','Latex', ...
                    'FontSize',16 );
                xlabel(axs_set, '$$ x $$', 'Interpreter','Latex','FontSize',16);
                yname = @(i_,k_) ['$$ d_x^' num2str(k_) 'u_' num2str(i_)  '$$'];
                xname_3D = '$$ x $$';
                yname_3D = '$$ u $$';
                zname_3D = '$$ d_x u $$';
            end

            if (ndim==3)
                axik = @(i_,k_) axs(max([i_,k_]));
            else
                if ((size(axs_mat,1)==1)||(size(axs_mat,2)==1))
                    axik = @(i_,k_) axs(max([i_,k_]));
                else
                    axik = @(i_,k_) axs_mat(i_,k_);
                end
            end

            if (iscell(S_))
                Smat = ldaux.Scell_2_Smat(S_,ndim);
            elseif ( length(size(S_)) == 2 )
                Smat = S_;
            else
                Smat = reshape( S_, ndim, [] );
            end

            xvec = Smat(1,:);
            utns = reshape(Smat(2:end,:),ndep,eor+1,[]);

            ileg = 1;
            for k = 1:(eor+1)
                for i = 1:ndep
                    leg_out(ileg) = plot( axik(i,k), ...
                    xvec, reshape(utns(i,k,:),size(xvec)), ...
                    'LineStyle', 'none', ...
                    'LineWidth', 0.5, ...
                    'Color', pspc.Color, ...
                    'Marker', pspc.Marker, ...
                    'MarkerSize', pspc.MarkerSize, ...
                    'MarkerFaceColor', pspc.Color, ...
                    'MarkerEdgeColor', [0 0 0] ...
                    );
                    ylabel(axik(i,k), yname(i,k-1), ...
                     'Interpreter','Latex','FontSize',16);
                     ileg = ileg + 1;
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
                xlabel(axs(3), xname_3D, 'Interpreter','Latex','FontSize',16);
                ylabel(axs(3), yname_3D, 'Interpreter','Latex','FontSize',16);
                zlabel(axs(3), zname_3D, 'Interpreter','Latex','FontSize',16);
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
                                plot(axik(i,k), ...
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
            if (~isfield(spc_out,'MarkerSize'))
                spc_out.MarkerSize = 8;
            end
            if (~isfield(spc_out,'LineStyle'))
                spc_out.LineStyle = '-';
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
        function obj_out = set_screen_posdim(obj,grid_dim_,plot_dim_,origin_tile_,screen_)
            if (nargin==2)
                posdim_use = grid_dim_;
                obj_out = obj;
                obj_out.fig = figure( ...
                'Name',obj.name, ...
                'WindowStyle','normal', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'Theme', 'dark', ...
                'Position', grid_dim_ ...
                );
            else
                if (nargin==5)
                    grid_dim = grid_dim_;
                    plot_dim = plot_dim_;
                    origin_tile = origin_tile_;
                    screen = screen_;
                else (nargin == 4)
                    grid_dim = grid_dim_;
                    plot_dim = plot_dim_;
                    origin_tile = origin_tile_;
                    screen = 1;
                end

                sys_screens = apv_plots.get_sys_screens();

                if (screen>size(sys_screens,1))
                    screen_i = sys_screens(1,:); % default to screen 1
                else
                    screen_i = sys_screens(screen,:);
                end

                o_screen_i = screen_i(1:2); d_screen_i = screen_i(3:4)-1;
                dels_grid_i = (d_screen_i)./[grid_dim(2) grid_dim(1)];
                dels_plot_i = dels_grid_i.*[plot_dim(2) plot_dim(1)];

                plt_lpos = floor( o_screen_i(1) + dels_grid_i(1)*( origin_tile(2)-1 ) );
                plt_bpos = floor( o_screen_i(2) + dels_grid_i(2)*( grid_dim(1) - origin_tile(1) ) );
                plt_wlen = floor(dels_plot_i(1));
                plt_hlen = floor(dels_plot_i(2));

                pos_set = [plt_lpos plt_bpos plt_wlen plt_hlen];

                fig_out = figure( ...
                'Name',obj.name, ...
                'WindowStyle','normal', ...
                'MenuBar', 'none', ...
                'ToolBar', 'none', ...
                'Theme', 'dark', ...
                'Units', 'pixels', ...
                'AutoResizeChildren', 'on' ...
                );
                set(fig_out,'OuterPosition',pos_set);
                obj_out = obj;
                obj_out.fig = fig_out;
            end
        end

        function sys_screens_out = get_sys_screens()
            sys_screens_out = get(groot,'MonitorPositions');
            if (~ismac) % works fine if osx
                arch = getenv('ARCH');
                [istart,iend] = regexp(arch,'mac');
                if ( (length(istart)*length(iend)) == 0 ) % works fine if osx
                    [istart,iend] = regexp(arch,'win'); % assume works fine if windows
                    if ( (length(istart)*length(iend)) == 0 ) % assume linux
                        % ScreenPixelsPerInch = java.awt.Toolkit.getDefaultToolkit().getScreenResolution()
                        ScreenDevices = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getScreenDevices();
                        MainScreen = java.awt.GraphicsEnvironment.getLocalGraphicsEnvironment().getDefaultScreenDevice().getScreen()+1;
                        MainBounds = ScreenDevices(MainScreen).getDefaultConfiguration().getBounds();
                        MonitorPositions = zeros(numel(ScreenDevices),4);
                        for n = 1:numel(ScreenDevices)
                            Bounds = ScreenDevices(n).getDefaultConfiguration().getBounds();
                            MonitorPositions(n,:) = [Bounds.getLocation().getX() + 1,-Bounds.getLocation().getY() + 1 - Bounds.getHeight() + MainBounds.getHeight(),Bounds.getWidth(),Bounds.getHeight()];
                        end
                        sys_screens_out = MonitorPositions;
                    end
                end
            end
        end
    end
end
