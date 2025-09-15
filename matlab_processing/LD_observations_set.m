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
        sigma_s;
        sigma_dnp1xu;
        sigma_JFs;

        crvs;
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
            elseif (nargin == 6)
                obj.eqn_name = eqn_;
                obj.nse_name = nse_;
                obj.gen_name = gen_;
                % dat_name = asmbl_gen_name(eqn_,nse_,gen_);
                % dat_suff = fam_;
                dat_name_adtl = fam_;
                dat_name = [asmbl_gen_name(eqn_,nse_,gen_) dat_name_adtl];
                dat_suff = bor_;
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

            %% conditional loading of extra data
            obj = obj.load_JFs;
            obj = obj.load_dnp1xu;
            obj = obj.load_sigmas;

            obj.crvs = [];
        end
        function cell_out = read_Sobs_cell(obj,xtra_)
            pts_struct = read_pts_struct([obj.dir_name '/' obj.dat_name xtra_ '.' obj.dat_suff]);
            cell_out = LD_observations_set.dat_cell(obj.ndim, ...
                                                        pts_struct.pts_in, ...
                                                        1:(pts_struct.ncrv), ...
                                                        pts_struct.npts_per_crv);
        end
        function obj_out = read_Sobs_slice(obj,xtra_)
            obj_out = obj;
            obj_out.reconstructed_set = true;
            dat_name = [obj.dat_name, xtra_];
            name_ = [obj.dir_name '/' dat_name '.' obj.dat_suff];

            pts_struct = read_pts_struct(name_);

            obj_out.dat_name = dat_name;

            obj_out.name = name_;

            % obj_out.eor = pts_struct.eor;
            % obj_out.ndep = pts_struct.ndep;
            % obj_out.ncrv = pts_struct.ncrv;
            % obj_out.nobs = pts_struct.nobs;
            obj_out.npts_per_crv = pts_struct.npts_per_crv;
            obj_out.pts_in = pts_struct.pts_in;
            % obj_out.ndim = 1 + obj.ndep*(obj.eor+1);
            obj_out.npts_uniform = prod(double(obj_out.npts_per_crv==obj_out.npts_per_crv(1)));

            obj_out.pts_mat = reshape(obj_out.pts_in,obj_out.ndim,[]);
        end
        function dnse_summary = read_denoise_summary(obj,sum_name_)

            name = [obj.dir_name '/' obj.dat_name sum_name_ '.' obj.dat_suff];
            file = fopen(name);
            if (file == -1)
                fprintf('(read_jet_sol_h_data) : WARNING - failed to read %s. Defaulting to empty output \n',name);
                dnse_summary = [];
            else
                hlen = fread(file,1,'int=>int');
                header = fread(file,hlen,'int=>int');
                [nwrite,nsmooth,ilen,dlen] = deal(header(1),header(2),header(3),header(4));
                idata = fread(file,ilen,'int=>int');
                ddata = fread(file,dlen,'double=>double');
                fclose(file);

                iwrite = idata(1:nwrite);
                ranks = idata((nwrite+1):end);
                residuals = ddata(1:nsmooth);

                dnse_summary = struct(  'nwrite', nwrite, ...
                                        'nsmooth', nsmooth, ...
                                        'iwrite', iwrite, ...
                                        'ranks', ranks, ...
                                        'residuals', residuals );
            end
        end
        function dxuk_out = read_dxuk_data(obj,name_)
            dxuk_out = read_dxuk_struct([obj.dir_name '/' obj.dat_name name_ '.' obj.dat_suff]);
        end
        function [R_svd,R_h_svd,pSj] = read_jet_sol_h_data(obj,nSVD_,nSmat_,ps_)
            if (nargin == 4)
                ps = ps_;
            else
                ps = '';
            end

            R_svd = obj.read_LD_svd([nSVD_{1} ps]);
            R_h_svd = obj.read_LD_svd([nSVD_{2} ps]);

            % theta_mat = LD_aux.read_Tmat([obj.dir_name '/' obj.dat_name nt_ '.' obj.dat_suff]);

            pSj = read_pts_struct([obj.dir_name '/' obj.dat_name nSmat_{1} ps '.' obj.dat_suff]);
            for i = 2:length(nSmat_)
                pSj(i) = read_pts_struct([obj.dir_name '/' obj.dat_name nSmat_{i} ps '.' obj.dat_suff]);;
            end

            for i = 1:length(nSmat_)
                pSj(i).pts_crv_inds = LD_observations_set.pts_crv_inds(obj.ndim,pSj(i).npts_per_crv);
            end
        end
        function rowimg_out = read_rowspace_image(obj,name_,fam_,bor_)
            name_tc = [obj.dir_name '/' obj.dat_name '_' fam_ '.' num2str(bor_) '.' name_ '_theta_chunk.' obj.dat_suff];
            theta_mat = read_Tmatrix_file(name_tc, 'double');

            name_vc = [obj.dir_name '/' obj.dat_name '_' fam_ '.' num2str(bor_) '.' name_ '_vnu_chunk.' obj.dat_suff];
            vnu_mat = read_Tmatrix_file(name_vc, 'double');

            name_vsc = [obj.dir_name '/' obj.dat_name '_' fam_ '.' num2str(bor_) '.' name_ '_vnu_syn_chunk.' obj.dat_suff];
            vnu_syn_mat = read_Tmatrix_file(name_vsc, 'double');

            rowimg_out = struct('thlmat', theta_mat', ...
                                'vnumat', vnu_mat', ...
                                'vnusmat', vnu_syn_mat');
        end
        function svd_out = read_LD_svd(obj,name_,fam_,bor_)
            switch nargin
                case 2
                    name = [obj.dir_name '/' obj.dat_name '.' name_ '.' obj.dat_suff];
                case 4
                    name = [obj.dir_name '/' obj.dat_name '_' fam_ '.' num2str(bor_) '.' name_ '.' obj.dat_suff];
            end
            svd_cell = read_Tmatrix_file(name,'double');
            if (length(svd_cell) == 3)
                [U,s,V] = deal(svd_cell{1},svd_cell{2},svd_cell{3});
            else
                [U,s,V] = deal([],svd_cell{1},svd_cell{2});
            end
            svd_out = struct('U',U,'s',s,'V',V);
        end
        function obj_out = load_sigmas(obj,suffix_)
            if (nargin==1)
                suffix = 'sigmas';
            else
                suffix = suffix_;
            end
            obj_out = obj;
            name = [obj.dir_name '/' obj.dat_name '_' suffix '.' obj.dat_suff];
            file = fopen(name);
            if (file == -1)
                [sigma_s,sigma_dnp1xu,sigma_JFs] = deal([]);
            else
                hlen = fread(file,1,'int=>int');
                header = fread(file,hlen,'int=>int');
                [eor,ndep] = deal(header(1),header(2));
                sigma_s = fread(file,1 + ndep*(eor+1),'double=>double');
                sigma_dnp1xu = fread(file,ndep,'double=>double');
                sigma_JFs = fread(file,ndep*(1 + ndep*(eor+1)),'double=>double');
                fclose(file);
            end
            obj_out.sigma_s = sigma_s;
            obj_out.sigma_dnp1xu = sigma_dnp1xu;
            obj_out.sigma_JFs = sigma_JFs;
        end
        function obj_out = load_dnp1xu(obj,suffix_)
            if (nargin==1)
                suffix = 'dnp1xu';
            else
                suffix = suffix_;
            end
            obj_out = obj;
            obj_out.dnp1xu_in = obj.read_adtl(suffix);
        end
        function obj_out = load_JFs(obj,suffix_)
            if (nargin==1)
                suffix = 'JFs';
            else
                suffix = suffix_;
            end
            obj_out = obj;
            obj_out.JFs_in = obj.read_adtl(suffix);
        end
        function adtl_out = read_adtl(obj,suffix_)
            suffix = suffix_;
            name = [obj.dir_name '/' obj.dat_name '_' suffix '.' obj.dat_suff];
            file = fopen(name);
            if (file == -1)
                fprintf('(read_adtl) : WARNING - failed to read %s. Defaulting to empty output \n',name);
                adtl_out = [];
            else
                hlen = fread(file,1,'int=>int');
                header = fread(file,hlen,'int=>int');
                adtl_in = fread(file,header(1),'double=>double');
                adtl_out = adtl_in;
                fclose(file);
            end
        end
        function [mats_out,svd_cell] = read_encoded_matrices(obj,mname_)
            mname = mname_;

            name = [obj.dir_name '/' obj.dat_name '_' mname '.' obj.dat_suff];
            file = fopen(name);
            hlen = fread(file,1,'int=>int');
            header = fread(file,hlen,'int=>int');
            [nset,ncol_full,ncon] = deal(header(1),header(2),header(3));
            mats_out = LD_encoded_matrix.empty(nset,0);

            for i = 1:nset
                dims_i = fread(file,2,'int=>int');
                vec_i = fread(file,prod(dims_i),'double=>double');
                mats_out(i) = LD_encoded_matrix(dims_i,reshape(vec_i,dims_i(2),dims_i(1))');
            end
            fclose(file);

            if (nargout == 2)
                name = [obj.dir_name '/' obj.dat_name '_' mname 'svd' '.' obj.dat_suff];
                file_svd = fopen(name);
                hlen = fread(file_svd,1,'int=>int');
                header = fread(file_svd,hlen,'int=>int');
                [nspc,vlen_full,ilen,dlen] = deal(header(1),header(2),header(3),header(4));
                nV_spcvec = fread(file,nspc,'int=>int');
                if ((ilen-(2*nspc))~=(vlen_full*nspc))
                    iV_spcmat = nan(nspc,vlen_full);
                    for i = 1:nspc
                        iV_spcmat(i,1:nV_spcvec(i)) = fread(file_svd,nV_spcvec(i),'int=>int');
                    end
                else
                    iV_spcmat = double(0:(vlen_full-1)) .* ones(nspc,vlen_full);
                end
                rank_vec = fread(file_svd,nspc,'int=>int');
                svd_cell = cell([5,nspc]); % nV, iVvec, rank, Vmat, svec
                for i = 1:nspc
                    svd_cell{1,i} = nV_spcvec(i);
                    svd_cell{2,i} = iV_spcmat(i,1:nV_spcvec(i));
                    svd_cell{3,i} = rank_vec(i);
                    svd_cell{4,i} = fread(file_svd,[nV_spcvec(i) nV_spcvec(i)],'double=>double');
                end
                for i = 1:nspc
                    svd_cell{5,i} = fread(file_svd,nV_spcvec(i),'double=>double');
                end
                fclose(file_svd);
            end
        end
        function [curves_out,icrvs] = make_curve_array(obj,icrvs_)
            if (nargin == 2)
                icrvs = icrvs_;
            else
                icrvs = 1:(obj.ncrv);
            end
            len_icrvs = length(icrvs(:));
            curves_out = LD_curve.empty(len_icrvs,0);

            [eor,ndep,ndim] = deal(obj.eor,obj.ndep,obj.ndim);

            [dnp1xu_loaded,JFs_loaded] = deal(~isempty(obj.dnp1xu_in),~isempty(obj.JFs_in));
            load_flag = double(dnp1xu_loaded)+double(JFs_loaded);

            npts_per_crv = obj.npts_per_crv;
            pts_in = obj.pts_in;
            inds_pts = LD_observations_set.dat_crv_inds(ndim,npts_per_crv);
            if (logical(load_flag))
                [JFs_in,dnp1xu_in] = deal(obj.JFs_in,obj.dnp1xu_in);
                inds_dnp1xu = LD_observations_set.dat_crv_inds(ndep,npts_per_crv);
                inds_JFs = LD_observations_set.dat_crv_inds(ndep*ndim,npts_per_crv);
                if (load_flag==2) % n+1 derivatives AND JFs loaded
                    for i = 1:len_icrvs
                        icrv = icrvs(i);
                        curves_out(i) = LD_curve(   eor, ndep, ...
                            pts_in(inds_pts(1,icrv):inds_pts(2,icrv)), ...
                            dnp1xu_in(inds_dnp1xu(1,icrv):inds_dnp1xu(2,icrv)), ...
                            JFs_in(inds_JFs(1,icrv):inds_JFs(2,icrv)) );
                    end
                elseif (dnp1xu_loaded) % just n+1 derivatives loaded
                    for i = 1:len_icrvs
                        icrv = icrvs(i);
                        curves_out(i) = LD_curve(   eor, ndep, ...
                            pts_in(inds_pts(1,icrv):inds_pts(2,icrv)), ...
                            dnp1xu_in(inds_dnp1xu(1,icrv):inds_dnp1xu(2,icrv)), ...
                            [] );
                    end
                else % just JFs loaded
                    for i = 1:len_icrvs
                        icrv = icrvs(i);
                        curves_out(i) = LD_curve(   eor, ndep, ...
                            pts_in(inds_pts(1,icrv):inds_pts(2,icrv)), ...
                            [], ...
                            JFs_in(inds_JFs(1,icrv):inds_JFs(2,icrv)) );
                    end
                end
            else % niether n+1 derivatives nor JFs loaded
                for i = 1:len_icrvs
                    icrv = icrvs(i);
                    curves_out(i) = LD_curve(eor,ndep, ...
                            pts_in(inds_pts(1,icrv):inds_pts(2,icrv)));
                end
            end

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
        function cell_out = dnp1xu_cell(obj,icrvs_)
            if (nargin == 2)
                icrvs = icrvs_;
            else
                icrvs = 1:(obj.ncrv);
            end
            cell_out = LD_observations_set.dat_cell(obj.ndep, ...
                                                    obj.dnp1xu_in, ...
                                                    icrvs, ...
                                                    obj.npts_per_crv);
        end
        function cell_out = JFs_cell(obj,icrvs_)
            if (nargin == 2)
                icrvs = icrvs_;
            else
                icrvs = 1:(obj.ncrv);
            end
            cell_out = LD_observations_set.dat_cell(obj.ndep*obj.ndim, ...
                                                    obj.JFs_in, ...
                                                    icrvs, ...
                                                    obj.npts_per_crv);
        end
        function cell_out = pts_cell(obj,icrvs_)
            if (nargin == 2)
                icrvs = icrvs_;
            else
                icrvs = 1:(obj.ncrv);
            end
            cell_out = LD_observations_set.dat_cell(obj.ndim, ...
                                                    obj.pts_in, ...
                                                    icrvs, ...
                                                    obj.npts_per_crv);
        end
        function pts_mat_out = pts_mat_crvi(obj,i_)
            [ncrv,ndim,pts_in] = deal(obj.ncrv,obj.ndim,obj.pts_in);
            inds = LD_observations_set.pts_crv_inds(ndim,obj.npts_per_crv);
            pts_mat_out = reshape(pts_in(inds(1,i_):inds(2,i_)),ndim,[]);
        end
        function dnp1xu_out = dnp1xu_mat(obj)
            dnp1xu_out = reshape(obj.dnp1xu_in,obj.ndep,[]);
        end
    end
    methods (Static)
        function pts_cell_out = combine_pts_cells(pSs01_,pSref_)
            if (nargin == 2) % lazy update
                pts_cell_out = pSref_;
                for i = 1:length(pSref_)
                    pts_cell_out{i}(:,1) = 0.5*(pSref_{i}(:,1) + pSs01_{i,1}(:,1));

                    pts_cell_out{i}(:,2:(end-1)) = (pSref_{i}(:,2:(end-1)) + ...
                                                    pSs01_{i,1}(:,2:end ) + ...
                                                    pSs01_{i,2}(:,1:(end-1) ) )/3.0;

                    pts_cell_out{i}(:,end) = 0.5*(pSref_{i}(:,end) + pSs01_{i,2}(:,end));
                end
            end
        end
        function pts_cell_out = pts_struct_2_cell(pS_,icrvs_)
            if (nargin == 2)
                icrvs = icrvs_;
            else
                icrvs = 1:(pS_(1).ncrv);
            end
            pS1 = pS_(1);

            ncell = prod(size(pS_));
            pts_cell_1 = LD_observations_set.dat_cell( 1+((pS1.ndep)*(pS1.eor+1)), ...
                                                    pS1.pts_in, ...
                                                    icrvs, ...
                                                    pS1.npts_per_crv  );
            if (ncell==1)
                pts_cell_out = pts_cell_1;
            else
                pts_cell_out = cell([pS_(1).ncrv,ncell]);
                pts_cell_out(:,1) = pts_cell_1;
                for i = 2:ncell
                    pts_cell_out(:,i) = LD_observations_set.dat_cell( 1+((pS_(i).ndep)*(pS_(i).eor+1)), ...
                                                            pS_(i).pts_in, ...
                                                            icrvs, ...
                                                            pS_(i).npts_per_crv  );
                end
            end

        end
        function crvs_out = regularize_curve_jets(crvs_,sigma_,lam_)
            if (nargin==1)
                sigma = ones(crvs_(1).jets.kor+1,crvs_(1).jets.ndep);
                lam = ones(crvs_(1).jets.ndep);
            else
                [sigma,lam] = deal(sigma_,lam_);
            end
            ncrvs = length(crvs_);
            crvs_out = crvs_;

            for i = 1:ncrvs
                crvs_out(i).rjets = LD_jets.regularize_jet(crvs_(i).jets,sigma,lam);
            end
        end
        function dat_cell_out = dat_cell(dlen_,dat_in_,icrvs_,npts_per_crv_)
            icrvs = icrvs_;
            len_icrvs = length(icrvs(:));
            inds = LD_observations_set.dat_crv_inds(dlen_,npts_per_crv_);
            dat_cell_out = cell([len_icrvs,1]);
            for i_cell = 1:len_icrvs
                icrv = icrvs(i_cell);
                dat_cell_out{i_cell,1} = reshape(dat_in_(inds(1,icrv):inds(2,icrv)), ...
                dlen_,[]);
            end
        end
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
            inds_out = LD_observations_set.dat_crv_inds(ndim_,npts_per_crv_);
        end
        function inds_out = dat_crv_inds(dlen_,npts_per_crv_)
            ncrv = length(npts_per_crv_);
            chunk_len_vec = dlen_*npts_per_crv_;
            inds_out = nan(2,ncrv);
            i_start = 0;
            for i = 1:ncrv
                inds_out(:,i) = [ i_start + 1; i_start + chunk_len_vec(i) ];
                i_start = i_start + chunk_len_vec(i);
            end
        end
    end
end
function mat_data_out = read_Tmatrix_file(name_,Tstr_)
    if (nargin == 1)
        Tstr = 'double';
    else
        Tstr = Tstr_;
    end
    T2T = [Tstr '=>' Tstr];
    file = fopen(name_);
    if (file == -1)
        fprintf('(read_Tmatrix_file) : ERROR - failed to read %s \n',name_);
        cell_out = 0;
    else
        hlen = fread(file,1,'int=>int');
        header = fread(file,hlen,'int=>int');
        switch hlen
            case 1 % implies multiple matrices
                nmat = header(1);
                mat_data_out = cell([nmat,1]);
                for i = 1:nmat
                    head_i = fread(file,2,'int=>int');
                    mat_data_out{i} = (fread(file,[head_i(2) head_i(1)],T2T))';
                end
            case 2 % implies just one matrix
                mat_data_out = (fread(file,[header(2) header(1)],T2T))';
        end
        fprintf('(read_Tmatrix_file) : read %s \n',name_);
        fclose(file);
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
function dxuk_struct = read_dxuk_struct(name_)
    file = fopen(name_);
    if (file == -1)
        fprintf('(read_dxuk_struct) : ERROR - failed to read %s \n',name_);
        dxuk_struct = 0;
    else
        ode_meta = fread(file,LD_observations_set.ode_meta_len,'int=>int');
        [kor,ndep] = deal(ode_meta(1),ode_meta(2));
        obs_meta = fread(file,LD_observations_set.obs_meta_len,'int=>int');
        [ncrv,nobs] = deal(obs_meta(1),obs_meta(2));
        npts_per_crv = fread(file,ncrv,'int=>int');
        dxuk_in = fread(file,kor*nobs,'double=>double');

        dxuk_struct = struct(    'kor', kor, ...
                                'ndep', ndep, ...
                                'ncrv', ncrv, ...
                                'nobs', nobs, ...
                                'npts_per_crv', npts_per_crv, ...
                                'dxuk_in', dxuk_in);
        fclose(file);
        fprintf('(read_dxuk_struct) : read %s \n',name_);
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

        pts_struct = struct(    'eor', eor, ...
                                'ndep', ndep, ...
                                'ncrv', ncrv, ...
                                'nobs', nobs, ...
                                'npts_per_crv', npts_per_crv, ...
                                'pts_in', pts_in);
        fclose(file);
        fprintf('(read_pts_struct) : read %s \n',name_);
    end
end
