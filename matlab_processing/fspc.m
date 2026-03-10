classdef fspc < jspc

    properties

        bor;

    end

    methods
        function obj = fspc(obs_)

            obj@jspc(obs_.ndep,obs_.eor);

        end
    end

end
