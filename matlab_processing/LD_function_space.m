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
            obj.bor = double(bor_);
            obj.perm_len = double(perm_len_);

            obj.eor = double(meta_.eor);
            obj.ndep = double(meta_.ndep);
            obj.ndim = double(meta_.ndim);

            obj.ndof_full = (1 + obj.ndep)*obj.perm_len;
        end
    end
end
