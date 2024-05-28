classdef LD_aux
    properties

    end
    methods
        function obj = LD_aux()

        end
    end
    methods (Static)
        function Asvd_rstr_out = comp_restricted_space(meta_,Lsvd_,Asvd_)
            Asvd_rstr_out = Asvd_;

            nvar = meta_.ndep + 1;
            ncrv = length(Lsvd_.rvec);
            ndof = Asvd_.ncol;
            rho_L = min(Lsvd_.rvec);
            kappa_A = ndof - max(Asvd_.rvec);
            ncols_PL = nvar*rho_L;

            Atns = permute(reshape(Asvd_.matT,ndof,[],ncrv),[2 1 3]);

            AYtns = LD_aux.Atns_Ytns_mult(Atns,Lsvd_.Vtns(:,1:rho_L,:));
            VYtns = LD_aux.Atns_Ytns_mult(Asvd_.Vtns,Lsvd_.Vtns(:,1:rho_L,:));

            rvec_AY = nan(ncrv,1);
            smat_AY = nan(ncols_PL,ncrv);
            Vtns_AY = nan(ncols_PL,ncols_PL,ncrv);
            for i = 1:ncrv
                rvec_AY(i) = rank(AYtns(:,:,i));
                [~,smat_AY(:,i),Vtns_AY(:,:,i)] = svd(AYtns(:,:,i),'econ','vector');
            end

            Asvd_rstr_out.AYtns = AYtns;
            Asvd_rstr_out.rvec_AY = rvec_AY;
            Asvd_rstr_out.smat_AY = smat_AY;
            Asvd_rstr_out.Vtns_AY = Vtns_AY;
            Asvd_rstr_out.VYtns = VYtns;
        end
        function AYtns_out = Atns_Ytns_mult(Atns_,Ytns_)
            [nrow,ndof,ncrv] = size(Atns_);
            [perm_len,rho,~] = size(Ytns_);
            nvar = ndof/perm_len;
            N_PL = nvar*rho;

            AYtns_out = nan(nrow,N_PL,ncrv);
            for i = 1:ncrv
                AYtns_out(:,:,i) = reshape(pagemtimes(reshape(Atns_(:,:,i),nrow,perm_len,nvar),Ytns_(:,:,i)),nrow,[]);
            end
        end
        function str1_out = overwrite_struct(str1_,str2_)
            str1_out = str1_;
            str1_names = fieldnames(str1_);
            for i = 1:length(str1_names)
                name_i = str1_names{i};
                str1_out = setfield(str1_out,name_i,getfield(str2_,name_i));
            end
        end
        function struct_out = combine_structs(str1_,str2_)
            struct_out = cell2struct([struct2cell(str1_); struct2cell(str2_)], [fieldnames(str1_); fieldnames(str2_)]);
        end
    end
end
