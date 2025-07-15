classdef LD_encoded_matrix
    properties (Constant)

    end
    properties

        dims;
        mat;

    end
    methods
        function obj = LD_encoded_matrix(dims_,mat_)
            dims = dims_;
            mat = mat_;

            obj.dims = dims;
            obj.mat = mat;
        end
    end
end
