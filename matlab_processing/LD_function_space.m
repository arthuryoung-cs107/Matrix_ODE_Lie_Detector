classdef LD_function_space
    properties
        bor;
        perm_len;

        eor;
        ndep;
        ndim;

        ndof_full;
    end
    methods
        function obj = LD_function_space(bor_,perm_len_,meta_)
            obj.bor = bor_;
            obj.perm_len = perm_len_;

            obj.eor = meta_.eor;
            obj.ndep = meta_.ndep;
            obj.ndim = meta_.ndim;

            obj.ndof_full = (1 + obj.ndep)*perm_len_;
        end
    end
end
