classdef adlam < adobj

    properties (Constant)

        A = @(o_) o_.Jac ;
        Omap = @(o_) o_.Jac * o_.xu + o_.b;
        vevl = @(o_,v_) o_.Jac * v_(:);

    end

    properties

        xu;
        b;

    end

    methods
        function obj = adlam(xu_,A_,b_)
            obj@adobj( A_*(xu_(:)) + b_(:) , reshape(A_,length(b_(:)),length(xu_(:))) );
            obj.xu = xu_(:);
            obj.b =  b_(:);
        end

    end

    methods (Static)
        function pr1l_out = pr1(l_)
            ndep = length(l_.xu)-1;

            vo =
            Jo =

            pr1l_out = adobj(vo,Jo);

        end

        function vec_out = coordgrads2Jac(mat_)
            [M,nobs] = size(mat_);
            N = size( mat_(1).Jac,2 );
            for iobs = 1:nobs
                val_i = nan(M,1);
                Jac_i = nan(M,N);
                for i = 1:M
                    val_i(i) = mat_(i,iobs).val;
                    Jac_i(i,:) = mat_(i,iobs).Jac;
                end
                vec_out(iobs) = adobj(val_i,Jac_i);
            end
        end
        function obj_o = trunc_to_adobj(obj)
            obj_o = adobj(obj.val,obj.Jac);
        end
    end

end
