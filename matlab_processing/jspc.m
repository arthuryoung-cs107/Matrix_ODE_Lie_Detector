classdef jspc

    properties

        ndep;
        eor;

    end

    methods

        function obj = jspc(ndep_,eor_)

            obj.ndep = ndep_;
            obj.eor = eor_;

        end

    end

    methods (Static)

        function [Smat,nobs,nset,kor_obs,ndim_obs] = unpack_Scell(Scell_,ndep_)

            nset = prod(size(Scell_));
            ndim_obs = size(Scell_{1},1);
            kor_obs = (ndim_obs-1)/ndep_ - 1;
            Smat = jspc.Scell_2_Smat( Scell_,ndim_obs );
            nobs = size(Smat,2);

        end
        function Smat_out = Scell_2_Smat(Scell_,ndim_)
            if (nargin==1)
                ndim = size(Scell_{1},2);
            else
                ndim = ndim_;
            end
            Smat_out = reshape(cell2mat(cellfun( @(s_) s_(:) ,Scell_,'UniformOutput',false )),ndim,[]);
        end
        function [xu_out,nobs_out] = xumat_nobs(obj,xu_)
            xu_out = reshape(xu_,jspc.nvar(obj),[]);
            nobs_out = size(xu_out,2);
        end
        function ndim_out = ndim(obj)
            ndim_out = 1+obj.ndep*(obj.eor+1);
        end
        function nvar_out = nvar(obj)
            nvar_out = 1+obj.ndep;
        end
        function [ndep_,eor_,ndim_,nvar_] = unpack_jetspace_dims(obj)
            ndep_ = obj.ndep;
            eor_ = obj.eor;
            ndim_ = 1 + ndep_*(eor_+1);
            nvar_ = 1 + ndep_;
        end
    end

end
