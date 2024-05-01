classdef LD_power_space < LD_function_space
    properties
        order_mat;
        dim0_lens;
        not_D_zero_inds;
    end
    methods
        function obj = LD_power_space(bor_,meta_)
            [perm_len order_mat dim0_lens] = count_set_perm_len(bor_,meta_.ndep);
            obj@LD_function_space(bor_,perm_len,meta_);

            obj.order_mat = order_mat;
            obj.dim0_lens = dim0_lens;
            obj.not_D_zero_inds = ~(obj.order_mat==0);
        end
    end
end

function [perm_len_out pow_mat_out dim0_lens_out] = count_set_perm_len(bor_,ndep_)
    d_ = ndep_ + 1;
    perm_len_out = (double(bor_+1))^(double(d_));

    dim0_lens_out = nan(bor_+1,1);
    pow_mat_out = nan(perm_len_out,d_);

    k_perm = 0;
    for i = 0:bor_
        [delk pow_mat_out] = set_powmat_full_recursive(pow_mat_out,bor_,d_,1,k_perm);
        dim0_lens_out(i+1) = delk;
        for k = k_perm:(k_perm+delk-1)
            pow_mat_out(k+1,1) = i;
        end
        k_perm = k_perm + delk;
    end
end

function [delk_out pow_mat_out] = set_powmat_full_recursive(mat_, bor_, d_, ilevel_, k_perm_)
    pow_mat_out = mat_;
    if (ilevel_>=d_)
        delk_out = 0;
    else
        delk_level = 0;
        if (ilevel_ == (d_-1))
            for i = 0:bor_
                pow_mat_out(k_perm_+delk_level+1,ilevel_+1) = i;
                delk_level = delk_level + 1;
            end
        else
            for i = 0:bor_
                [delk pow_mat_out] = set_powmat_full_recursive(pow_mat_out,bor_,d_,ilevel_+1,delk_level + k_perm_);
                for k = 0:(delk-1)
                    pow_mat_out(k_perm_+delk_level+1,ilevel_+1) = i;
                    delk_level = delk_level + 1;
                end
            end
        end
        delk_out = delk_level;
    end
end
