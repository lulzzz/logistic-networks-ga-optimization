classdef NetworkSimulator
    properties
        simTime;
        n;
        m;
        L;
        LT;
        LA;
        LA_nom;
        d;
        B_nom;
        
        initialHC;
        alpha;
        beta;

        u; % init
        u_hist; % init
        HC = 0; % init
        h; % init
        x_inf = 100000; % infinite stock at external sources

        xd; % key parameter
        x; % key parameter

        isCrashed = false;
        satisfiedRate = 0.0;
        fitness = 0.0;
    end
    methods
        function obj = NetworkSimulator(simTime, n, m, L, LT, LA, d, LA_nom, B_nom)
            obj.simTime = simTime;
            obj.n = n;
            obj.m = m;
            obj.L = L;
            obj.LT = LT;
            obj.LA = LA;
            obj.LA_nom = LA_nom;
            obj.d = d;
            obj.B_nom = B_nom;

            obj.u = zeros(n, simTime+1);
            obj.u_hist = zeros(n, simTime+1);
        end

        function obj = simulate(obj, xd_min)
            obj.xd = [xd_min; 0; 0];
            obj.x = zeros(obj.n+obj.m, obj.simTime+1);
            xd_min_source(1:obj.m) = obj.x_inf;
            obj.x(1:obj.n+obj.m,1) = [xd_min; xd_min_source'];
            obj.isCrashed = false;
            obj.HC = 0;

            for t = 1:obj.simTime
                o_q = zeros(obj.n,1);
                for i=1:obj.L
                    for j=i:obj.L
                        if (t > obj.L) || (t - i > 0)
                            o_q = o_q + B(:,:,j,t - i)*obj.u_hist(:,t - i);
                        end
                    end
                end

                obj.u(:,t) = obj.xd(1:obj.n) - obj.x(1:obj.n,t) - o_q;

                obj.u(:,t) = max( obj.u(:,t), 0 );

                obj.u_hist(:,t) = obj.u(:,t);

                ur = zeros(obj.n,1);
                for j=1:obj.n
                    ur(j) = 0;
                    for i=1:obj.n+obj.m
                        if (t > obj.L) || (t - obj.LT(i,j) > 0)
                            ur(j) = ur(j)+obj.LA(i,j,t - obj.LT(i,j))*obj.u_hist(j,t - obj.LT(i,j));
                        end
                    end
                end

                obj.h(:,t) = min ( obj.d(:,t), obj.x(1:obj.n,t) + ur );

                y(:,t) = obj.x(:,t);
                y(1:obj.n,t) = obj.x(1:obj.n,t) + ur - obj.h(:,t);

                S_u(:,t) = obj.LA_nom*obj.u(:,t);

                obj.LA(:,:,t) = obj.LA_nom;
                for i=1:obj.n
                    temp = y(i,t) - S_u(i,t);
                    if temp < 0
                        for j=1:obj.n
                            obj.LA(i,j,t) = obj.LA_nom(i,j)*(1 + temp/S_u(i,t));
                            if obj.LA(i,j,t) < 0
                                flaga = 1
                            end
                        end
                    end
                end

                % Update delay matrices B for delay 1 to L, update of B_0 is not required
                B(:,:,:,t) = obj.B_nom;
                for k = 1:obj.L % index k corresponds to delay k
                    for j = 1:obj.n
                        t_sum = 0;
                        for i = 1:obj.n+obj.m
                            if obj.LT(i,j) == k
                                t_sum = t_sum + obj.LA(i,j,t);
                            end
                        end
                        B(j,j,k,t) = t_sum;
                    end
                end
                S_u_mod(:,t) = obj.LA(:,:,t)*obj.u_hist(:,t);

                % Stock balance equation
                obj.x(:,t+1) = y(:,t) - S_u_mod(:,t);

                temp_d = obj.d(:,1:t);
                temp_h = obj.h(:,1:t);

                if sum(temp_h(:)) / sum(temp_d(:)) < 1
                    obj.isCrashed = true;
                end

                obj.HC = obj.HC + sum(obj.x(1:obj.n, t+1));
            end
            real_demand = obj.d(:,1:obj.simTime);
            obj.satisfiedRate = sum(obj.h(:)) / sum(real_demand(:));
            
            if obj.initialHC
                obj.fitness = (1 - (obj.HC / obj.initialHC))^obj.alpha * obj.satisfiedRate^obj.beta + 1;
            end
%             obj.fitness = 1+(10000000/(round(obj.HC)+1)) * power((1/(2 - obj.satisfiedRate)), 60);
        end
    end
end
