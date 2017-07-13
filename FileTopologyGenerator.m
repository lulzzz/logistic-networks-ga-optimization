classdef FileTopologyGenerator
    properties
        simTime;
        m;
        n;
        
        LT;
        L; % max lead time
        
        LA_nom;
        LA;
        
        d;
    end
    methods
        function obj = FileTopologyGenerator(filepath)
            topologyFile = load(filepath);
            obj.simTime = topologyFile.network.simTime;
            obj.m = topologyFile.network.m; % nr of external sources
            obj.n = topologyFile.network.n;  % nr of controlled nodes
            

            obj.LT = topologyFile.network.LT;            
            obj.L = max(obj.LT(:));
            obj.LA_nom(:,:,1) = topologyFile.network.LA_nom;
            obj.LA = zeros(obj.n+obj.m,obj.n,obj.simTime+1);
            obj.LA(:,:,1) = obj.LA_nom;
            
            % Verify if allocation correct - elements in each column should sum up to 1 or 0
            for j=1:obj.n
                temp = 0;
                for i=1:obj.n+obj.m
                    temp = temp + obj.LA(i,j,1);
                end
                if (temp == 0) || (temp == 1)
                    fprintf('Proper allocation in column: %d\n', j);
                else
                    error('Improper allocation in column: %d', j);
                end
            end
            
            obj.d = topologyFile.network.d;
        end
    end
end
