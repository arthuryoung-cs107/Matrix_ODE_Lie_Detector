classdef adfcn

    properties (Constant)

    end

    properties

        val; % scalar
        grad; % tall vector

    end

    methods
        function obj = adfcn(val_,grad_)
            [obj.val,obj.grad] = deal(val_(1),grad_(:));
        end

        %% overloaded operators

        function obj_o = plus(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn(vl+vr,gl+gr);
        end
        function obj_o = minus(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn(vl-vr,gl-gr);
        end
        function obj_o = uminus(obj_)
            [vi,gi] = adfcn.unpack_valgrad(obj_);
            obj_o = adfcn(-vi,-gi);
        end
        function obj_o = uplus(obj_)
            [vi,gi] = adfcn.unpack_valgrad(obj_);
            obj_o = adfcn(vi,gi);
        end
        function obj_o = times(obj_l_,obj_r_)
            [vl,gl,vr,gr] = adfcn.unpack_valgrad_pair(obj_l_,obj_r_);
            obj_o = adfcn( vl.*vr , gl.*vr + vl.*gr );
        end

    end

    methods (Static)

        function f_ad_out = f_ad( s_,f_ )
            s_ad = adfcn.seed_sol(s_);
            f_ad_out = f_(s_ad(:));
        end
        function objs_o = seed_sol(s_)
            ndim = length(s_(:));
            objs_o = adfcn.empty(ndim,0);
            Jac = eye(ndim);
            for i = 1:ndim
                objs_o(i) = adfcn( s_(i), Jac(:,i) );
            end
        end

        function [vl,gl,vr,gr] = unpack_valgrad_pair(obj_l_,obj_r_)
            if ( isnumeric(obj_l_) )
                [vl,gl] = deal(obj_l_,0);
            else
                [vl,gl] = deal(obj_l_.val,obj_l_.grad);
            end
            if ( isnumeric(obj_r_) )
                [vr,gr] = deal(obj_r_,0);
            else
                [vr,gr] = deal(obj_r_.val,obj_r_.grad);
            end
        end
        function [vo,go] = unpack_valgrad(obj_)
            if ( isnumeric(obj_) )
                [vo,go] = deal(obj_,0);
            else
                [vo,go] = deal(obj_.val,obj_.grad);
            end
        end
    end

end
