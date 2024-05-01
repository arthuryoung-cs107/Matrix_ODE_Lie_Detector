classdef LD_observations_set
    properties (Constant)
        ode_meta_len = 2;
        obs_meta_len = 2;
    end
    properties
        dir_name;
        dat_suff;
        dat_name;

        name;

        eor;
        ndep;
        ncrv;
        nobs;

        npts_per_crv;
        pts_in;

        ndim;

        JFs_in;
        dnp1xu_in;
    end
    methods
        function obj = LD_observations_set(dir_name_,dat_name_,dat_suff_)
            name_ = [dir_name_ '/' dat_name_ '.' dat_suff_];

            fprintf('(LD_observations_set::LD_observations_set) reading %s\n', name_);
            pts_struct = read_pts_struct(name_);

            obj.dir_name = dir_name_;
            obj.dat_name = dat_name_;
            obj.dat_suff = dat_suff_;

            obj.name = name_;

            obj.eor = pts_struct.eor;
            obj.ndep = pts_struct.ndep;
            obj.ncrv = pts_struct.ncrv;
            obj.nobs = pts_struct.nobs;
            obj.npts_per_crv = pts_struct.npts_per_crv;
            obj.pts_in = pts_struct.pts_in;

            obj.ndim = 1 + obj.ndep*(obj.eor+1);
        end
        function pts_cell_out = pts_cell(obj)
            [ncrv,ndim,pts_in] = deal(obj.ncrv,obj.ndim,obj.pts_in);
            inds = LD_observations_set.pts_crv_inds(ndim,obj.npts_per_crv);
            pts_cell_out = cell([ncrv,1]);
            for i = 1:ncrv
                pts_cell_out{i,1} = reshape(pts_in(inds(1,i):inds(2,i)),ndim,[]);
            end
        end
    end
    methods (Static)
        function meta_out = make_meta_data(eor_,ndep_)
            meta_out = struct(  'eor', eor_, ...
                                'ndep', ndep_, ...
                                'ndim', 1+ndep_*(eor_+1));
        end
        function inds_out = pts_crv_inds(ndim_,npts_per_crv_)
            ncrv = length(npts_per_crv_);
            chunk_len_vec = ndim_*npts_per_crv_;
            inds_out = nan(2,ncrv);
            i_start = 0;
            for i = 1:ncrv
                inds_out(:,i) = [ i_start + 1; i_start + chunk_len_vec(i) ];
                i_start = i_start + chunk_len_vec(i);
            end
        end
    end
end

function pts_struct = read_pts_struct(name_)
    file = fopen(name_);
    ode_meta = fread(file,LD_observations_set.ode_meta_len,'int=>int');
    [eor,ndep] = deal(ode_meta(1),ode_meta(2));
    obs_meta = fread(file,LD_observations_set.obs_meta_len,'int=>int');
    [ncrv,nobs] = deal(obs_meta(1),obs_meta(2));
    npts_per_crv = fread(file,ncrv,'int=>int');
    pts_in = fread(file,(1 + ndep*(eor+1))*nobs,'double=>double');
    fclose(file);

    pts_struct = struct(    'eor', eor, ...
                            'ndep', ndep, ...
                            'ncrv', ncrv, ...
                            'nobs', nobs, ...
                            'npts_per_crv', npts_per_crv, ...
                            'pts_in', pts_in);
end
