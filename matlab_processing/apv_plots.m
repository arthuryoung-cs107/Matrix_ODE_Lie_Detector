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
