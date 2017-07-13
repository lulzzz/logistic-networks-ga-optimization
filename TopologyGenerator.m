classdef TopologyGenerator
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
        function obj = TopologyGenerator(simTime, m)
            obj.simTime = simTime;
            obj.m = m; % nr of external sources
            obj.n = (obj.m*(obj.m+1))/2;  % nr of controlled nodes
            
            % Lead-time delays
            % LT(i,j) - delay from node i to j
            % Situation of supplying external source is excluded by restricting
            % dimension to (n+m,n)
            obj.LT = zeros(obj.n+obj.m,obj.n);
            
            for index=1:obj.n-1
               obj.LT(index, index+1) = round(1+10*rand(1,1));
            end

            for index=obj.n+1:obj.n+obj.m
               obj.LT(index, index-obj.n) = round(1+10*rand(1,1)); 
            end

            for index=obj.n:-1:obj.m+1
               obj.LT(index-2, index) = round(1+10*rand(1,1)); 
            end
            
            obj.L = max(obj.LT(:));
            
            % LA(i,j,k) - part of goods quantity requested from node i by j at instant k
            % Situation of supplying an external source is excluded by restricting the
            % dimension to (n+m,n)
            % LA_nom stores the nominal order partitioning coefficients
            obj.LA_nom = zeros(obj.n+obj.m,obj.n);
            obj.LA = zeros(obj.n+obj.m,obj.n,obj.simTime+1);
            for index=1:obj.n-1
               obj.LA_nom(index, index+1) = round(1+8*rand(1,1))/10;
            end

            for index=obj.n+1:obj.n+obj.m
               obj.LA_nom(index, index-obj.n) = 1 - sum(obj.LA_nom(:,index-obj.n));
            end

            for index=obj.n:-1:obj.m+1
               obj.LA_nom(index-2, index) = 1 - sum(obj.LA_nom(:,index));
            end
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
            
            % Random demand generator
            obj.d = zeros(obj.n,obj.simTime+1);
            
            for j=1:obj.simTime+1
                for node_index=1:obj.n
                    obj.d(node_index,j) = round( gamrnd(5,10) );
                end
            end
        end
    end
end
