classdef LD_observations_set
    properties (Constant)
        ode_meta_len = 2;
        obs_meta_len = 2;
    end
    properties
        reconstructed_set=false;
        eqn_name;
        nse_name;
        gen_name;
        fam_name;
        bor_name;
        mat_name;
        rec_name;

        nrow;
        ncol;
        matT;

        rvec;
        smat;
        Vtns;

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
        npts_uniform;

        pts_mat;

        JFs_in;
        dnp1xu_in;
    end
    methods
        function obj = LD_observations_set(dir_name_,eqn_,nse_,gen_,fam_,bor_,mat_,rec_,suf_,dat_suff_)
            if (nargin == 3)
                dat_name = eqn_;
                dat_suff = nse_;
            elseif (nargin == 5)
                obj.eqn_name = eqn_;
                obj.nse_name = nse_;
                obj.gen_name = gen_;
                dat_name = asmbl_gen_name(eqn_,nse_,gen_);
                dat_suff = fam_;
            else
                obj.reconstructed_set = true;
                obj.eqn_name = eqn_;
                obj.nse_name = nse_;
                obj.gen_name = gen_;
                obj.fam_name = fam_;
                obj.bor_name = bor_;
                obj.mat_name = mat_;
                obj.rec_name = rec_;
                [dat_name,mat_file_name] = asmbl_rec_name(eqn_,nse_,gen_,fam_,bor_,mat_,rec_,suf_);
                mat_svd_pckg = read_mat_svd_struct([dir_name_ '/' mat_file_name],dat_suff_);
                obj.nrow = mat_svd_pckg.nrow;
                obj.ncol = mat_svd_pckg.ncol;
                obj.matT = mat_svd_pckg.matT;
                obj.rvec = mat_svd_pckg.rvec;
                obj.smat = mat_svd_pckg.smat;
                obj.Vtns = mat_svd_pckg.Vtns;

                dat_suff = dat_suff_;
            end

            name_ = [dir_name_ '/' dat_name '.' dat_suff];

            fprintf('(LD_observations_set::LD_observations_set) reading %s\n', name_);
            pts_struct = read_pts_struct(name_);

            obj.dir_name = dir_name_;
            obj.dat_name = dat_name;
            obj.dat_suff = dat_suff;

            obj.name = name_;

            obj.eor = pts_struct.eor;
            obj.ndep = pts_struct.ndep;
            obj.ncrv = pts_struct.ncrv;
            obj.nobs = pts_struct.nobs;
            obj.npts_per_crv = pts_struct.npts_per_crv;
            obj.pts_in = pts_struct.pts_in;

            obj.ndim = 1 + obj.ndep*(obj.eor+1);
            obj.npts_uniform = prod(double(obj.npts_per_crv==obj.npts_per_crv(1)));

            obj.pts_mat = reshape(obj.pts_in,obj.ndim,[]);
        end
        function fspace_name_out = make_fspace_config_name(obj,fam_,bor_)
            fspace_name_out = [obj.dir_name '/' obj.dat_name '_' fam_ '.' num2str(bor_) '.domain_config.' obj.dat_suff];
        end
        function meta_out = meta_data(obj)
            meta_out = LD_observations_set.make_meta_data(obj.eor,obj.ndep);
        end
        function mat_pckg_out = read_mat_package(obj,fam_,bor_,mat_)
            mat_pckg_out = read_mat_struct([obj.dir_name '/' obj.dat_name '_' fam_ '.' num2str(bor_) '.' mat_],obj.dat_suff);
            mat_svd_pckg_out.dat_name = [obj.dat_name '_' fam_ '.' num2str(bor_) '.' mat_];
        end
        function [mat_svd_pckg_out,rstr_pckg_out] = read_mat_svd_package(obj,fam_,bor_,mat_,rstr_,mat_rstr_)
            base_name = [obj.dat_name '_' fam_ '.' num2str(bor_) '.'];
            if (nargin == 4)
                mat_svd_pckg_out = read_mat_svd_struct([obj.dir_name '/' base_name mat_],obj.dat_suff);
                mat_svd_pckg_out.dat_name = [base_name mat_];
            else
                Amat_pckg = read_mat_struct([obj.dir_name '/' base_name mat_],obj.dat_suff);
                Lmat_svd_pckg = read_mat_svd_struct([obj.dir_name '/' base_name mat_rstr_],obj.dat_suff);
                AYrstr_svd_pckg = read_mat_svd_struct([obj.dir_name '/' base_name mat_ rstr_],obj.dat_suff);
                if (nargout==2)
                    rstr_pckg_out = LD_observations_set.make_restricted_svd_package(Amat_pckg,Lmat_svd_pckg,AYrstr_svd_pckg);
                    mat_svd_pckg_out = AYrstr_svd_pckg;
                    mat_svd_pckg_out.nrow = rstr_pckg_out.nrow;
                    mat_svd_pckg_out.matT = rstr_pckg_out.matT;
                else
                    mat_svd_pckg_out = LD_observations_set.make_restricted_svd_package(Amat_pckg,Lmat_svd_pckg,AYrstr_svd_pckg);
                end
                mat_svd_pckg_out.dat_name = [base_name mat_ rstr_];
            end
        end
        function obj_out = read_additional_observations(obj,name1_)
            obj_out = obj;
            file = fopen(name1_);
            hlen = fread(file,1,'int=>int');
            header = fread(file,hlen,'int=>int');
            if (hlen==3)
                obj_out.dnp1xu_in = fread(file,header(1),'double=>double');
            elseif (hlen==4)
                obj_out.JFs_in = fread(file,header(1),'double=>double');
            end
            fclose(file);
        end
        function pckg_out = mat_package(obj)
            nconstr_dim = obj.nrow/obj.nobs;
            pckg_out = struct(  'mat', obj.matT', ...
                                'nconstr_dim', nconstr_dim, ...
                                'npts_per_crv', obj.npts_per_crv);
        end
        function pckg_out = mat_svd_pckg(obj)
            pckg_out = obj.mat_package;
            pckg_out.rvec = obj.rvec;
            pckg_out.smat = obj.smat;
            pckg_out.Vtns = obj.Vtns;
        end
        function pts_cell_out = pts_cell(obj)
            [ncrv,ndim,pts_in] = deal(obj.ncrv,obj.ndim,obj.pts_in);
            inds = LD_observations_set.pts_crv_inds(ndim,obj.npts_per_crv);
            pts_cell_out = cell([ncrv,1]);
            for i = 1:ncrv
                pts_cell_out{i,1} = reshape(pts_in(inds(1,i):inds(2,i)),ndim,[]);
            end
        end
        function pts_mat_out = pts_mat_crvi(obj,i_)
            [ncrv,ndim,pts_in] = deal(obj.ncrv,obj.ndim,obj.pts_in);
            inds = LD_observations_set.pts_crv_inds(ndim,obj.npts_per_crv);
            pts_mat_out = reshape(pts_in(inds(1,i_):inds(2,i_)),ndim,[]);
        end
    end
    methods (Static)
        function meta_out = make_meta_data(eor_,ndep_)
            meta_out = struct(  'eor', eor_, ...
                                'ndep', ndep_, ...
                                'ndim', 1+ndep_*(eor_+1));
        end
        function [AYLglbsvd_out,AYLglbsvd_glb_out] = make_global_restricted_svd_package(Asvd_,Lglbsvd_,kappa_L_)
            if (nargin == 3)
                kappa_L = kappa_L_;
            else
                kappa_L = 1;
            end

            AYLglbsvd_glb_out = LD_aux.make_global_restricted_svd_package(Asvd_,Lglbsvd_,kappa_L);
            AYLglbsvd_glb_out.dat_name = [Asvd_.dat_name 'YLglb(nofile)'];

            if (nargout == 2)
                AYLglbsvd_out = LD_aux.overwrite_struct(Asvd_,AYLglbsvd_glb_out);
            else
                AYLglbsvd_out = AYLglbsvd_glb_out;
            end
        end
        function [mat_svd_pckg_out, mat_rstr_svd_pckg_out] = make_restricted_svd_package(Amat_pckg_, Lmat_svd_pckg_, AYrstr_svd_pckg_)
            [ndof,ncrv,perm_len] = deal(size(Amat_pckg_.matT,1),size(Lmat_svd_pckg_.Vtns,3),size(Lmat_svd_pckg_.Vtns,1));
            nvar = ndof/perm_len;
            if (nargin==3)
                mat_rstr_svd_pckg_out = AYrstr_svd_pckg_;
                ncol = double(AYrstr_svd_pckg_.ncol);
                rho_PL = ncol/nvar;
                Ytns_L = Lmat_svd_pckg_.Vtns(:,1:rho_PL,:);
                AVLtns = LD_aux.Atns_Ytns_mult(permute(reshape(Amat_pckg_.matT,ndof,[],ncrv),[2 1 3]),Ytns_L);
            else
                mat_rstr_svd_pckg_out = Lmat_svd_pckg_;
                rho_PL = double(min(Lmat_svd_pckg_.rvec));
                ncol = rho_PL*nvar;
                Ytns_L = Lmat_svd_pckg_.Vtns(:,1:rho_PL,:);
                AVLtns = LD_aux.Atns_Ytns_mult(permute(reshape(Amat_pckg_.matT,ndof,[],ncrv),[2 1 3]),Ytns_L);
                rvec = nan(ncrv,1);
                smat = nan(ncol,ncrv);
                Vtns = nan(ncol,ncol,ncrv);
                for i = 1:ncrv
                    rvec(i) = rank(AVLtns(:,:,i));
                    [~,smat(:,i),Vtns(:,:,i)] = svd(AVLtns(:,:,i),'econ','vector');
                end
                mat_rstr_svd_pckg_out.ncol = ncol;
                mat_rstr_svd_pckg_out.rvec = rvec;
                mat_rstr_svd_pckg_out.smat = smat;
                mat_rstr_svd_pckg_out.Vtns = Vtns;
            end
            mat_rstr_svd_pckg_out.nrow = Amat_pckg_.nrow;
            mat_rstr_svd_pckg_out.matT = reshape(permute(AVLtns,[2 1 3]),ncol,[]);
            mat_rstr_svd_pckg_out.rho_PL = rho_PL;
            mat_rstr_svd_pckg_out.matT_raw = Amat_pckg_.matT;
            mat_rstr_svd_pckg_out.smat_L = Lmat_svd_pckg_.smat;
            mat_rstr_svd_pckg_out.Ytns_L = Ytns_L;
            if (nargout==2)
                mat_svd_pckg_out = LD_aux.overwrite_struct(Amat_pckg_,mat_rstr_svd_pckg_out);
                mat_svd_pckg_out.dat_name = [Amat_pckg_.dat_name 'YL(nofile)'];
            else
                mat_svd_pckg_out = mat_rstr_svd_pckg_out;
            end
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

function str_out = asmbl_gen_name(eqn_,nse_,gen_)
    str_out = [eqn_ '_' nse_ '_' gen_ 'gen'];
end
function [str_out,gen_name_out] = asmbl_mat_name(eqn_,nse_,gen_,fam_,bor_,mat_)
    gen_name_out = asmbl_gen_name(eqn_,nse_,gen_);
    str_out = [gen_name_out '_' fam_ '.' num2str(bor_) '.' mat_];
end
function [str_out,mat_name_out,gen_name_out] = asmbl_rec_name(eqn_,nse_,gen_,fam_,bor_,mat_,rec_,suf_)
    [mat_name_out,gen_name_out] = asmbl_mat_name(eqn_,nse_,gen_,fam_,bor_,mat_);
    str_out = [mat_name_out '.' rec_ suf_];
end
function mat_struct = read_mat_struct(name_,suf_)
    name_full = [name_ '.' suf_];
    file = fopen(name_full);
    if (file == -1)
        fprintf('(read_mat_struct) : ERROR - failed to read %s \n',name_full);
        mat_struct = 0;
    else
        hlen = fread(file,1,'int=>int');
        header = fread(file,hlen,'int=>int');
        if (hlen==2)
            [nrow,ncol] = deal(header(1),header(2));
            mat_in = fread(file,prod(header),'double=>double');
            mat_struct = struct(    'nrow', nrow, ...
                                    'ncol', ncol, ...
                                    'matT', reshape(mat_in,ncol,nrow));
        else
            fprintf('(read_mat_struct) : ERROR - attempted to read improperly formatted file from %s (hlen=%d) \n',name_full,hlen);
            mat_struct = 0;
        end
        fclose(file);
    end
end
function mat_svd_struct = read_mat_svd_struct(name_,suf_)
    name_full = [name_ '_svd.' suf_];
    file = fopen(name_full);
    if (file==-1)
        fprintf('(read_mat_struct) : ERROR - failed to read %s \n',name_full);
        mat_svd_struct = 0;
    else
        hlen = fread(file,1,'int=>int');
        header = fread(file,hlen,'int=>int');
        if (hlen==2)
            mat_svd_struct = read_mat_struct(name_,suf_);
            [ncrv,ncol] = deal(header(1),header(2));
            rvec = fread(file,ncrv,'int=>int');
            s_in = fread(file,ncrv*ncol,'double=>double');
            V_in = fread(file,ncrv*ncol*ncol,'double=>double');
            mat_svd_struct.rvec = rvec;
            mat_svd_struct.smat = reshape(s_in,ncol,ncrv);
            mat_svd_struct.Vtns = reshape(V_in,ncol,ncol,ncrv);
        elseif (hlen==3)
            [ncrv,ncol] = deal(header(1),header(3));
            rvec = fread(file,ncrv,'int=>int');
            s_in = fread(file,ncrv*ncol,'double=>double');
            V_in = fread(file,ncrv*ncol*ncol,'double=>double');
            mat_svd_struct = struct('nrow', [], ...
                                    'ncol', ncol, ...
                                    'matT', [], ...
                                    'rvec', rvec, ...
                                    'smat', reshape(s_in,ncol,ncrv), ...
                                    'Vtns', reshape(V_in,ncol,ncol,ncrv));
        else
            fprintf('(read_mat_svd_struct) : ERROR - attempted to read improperly formatted file from %s (hlen=%d) \n',name_full,hlen);
            mat_struct = 0;
        end
        fclose(file);
    end
end
function pts_struct = read_pts_struct(name_)
    file = fopen(name_);
    if (file == -1)
        fprintf('(read_pts_struct) : ERROR - failed to read %s \n',name_);
        pts_struct = 0;
    else
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
end
