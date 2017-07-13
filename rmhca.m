% logistic-networks-ga-optimization

clc;
clear all;

beginApplicationTic = tic();

chartVisibility = true;
fileGenerator = true;
generations = 5000;
stepRange = 10;
learningBreak = 500;

methodType = 'ga';
topologyFilename = 'topologies/2_100_a.mat';

% GA parameters
generationSize = 10;
mutationProbability = 0.01;
alpha = 10;
beta = 40;
resultGenerationNo = 0;

x_inf = 100000; % infinite stock at external sources
generationsDone = 0;

% TOPOLOGY GENERATOR
simTime = 1000;
m = 3;  % nr of external sources

% network = TopologyGenerator(simTime, m);
network = FileTopologyGenerator(topologyFilename);

simTime = network.simTime;
n = network.n;
m = network.m;
LT = network.LT;
L = network.L;
LA_nom = network.LA_nom;
LA = network.LA;
d = network.d;

% Initial conditions
time = linspace(0, simTime-1, simTime); % from 0 to simTime-1
u = zeros(n, simTime+1);
u_hist = zeros(n, simTime+1); % order history
x = zeros(n+m, simTime+1);
y = zeros(n+m, simTime+1);

dmax = max(d,[],2); % take the biggest value of each row

% State-space description
% System matrices
B_nom = zeros(n,n,L);
B = zeros(n,n,L,simTime+1);

% Assuming zero order processing time
B_0 = -LA(1:n,1:n,1);

for k = 1:L % index k corresponds to delay k
    for j = 1:n
        t_sum = 0;
        for i = 1:n+m
            if LT(i,j) == k
                t_sum = t_sum + LA(i,j,1);
            end
        end
        B_nom(j,j,k) = t_sum;
    end
end

B(:,:,:,1) = B_nom;

% Sum of delay matrices
Lambda = zeros(n);
for k = 1:L % table index k corresponds to delay k
    Lambda = Lambda + B(:,:,k,1);
end

Lambda = Lambda + B_0;

% Reference stock level for full demand satisfaction
temp = zeros(n);
for k = 1:L % table index k corresponds to delay k
    temp = temp + k*B(:,:,k,1);
end

dmean = mean(d,2);
dstd = std(d,0,2);
xd_min = ceil((eye(n) + temp)*inv(Lambda)*dmax) + 1

xd = [xd_min; 0; 0]; % reference stock level (see calculation of xd_min below)

xd_min_source(1:m) = x_inf;
x(1:n+m,1) = [xd_min; xd_min_source']; % initial stock level

start_xd_min = xd_min;
best_simulation = 0;
a = NetworkSimulator(simTime, n, m, L, LT, LA, d, LA_nom, B_nom);

% Main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = simulate(a, xd_min);
initial_HC = a.HC;

a.initialHC = initial_HC;
a.alpha = alpha;
a.beta = beta;

best_HC = a.HC;
best_xd_min = a.xd(1:n);
last_good_xdmin = best_xd_min;

methodTypeChars=char(methodType);
switch methodTypeChars
    case char('brute-force')
        for xd_min1 = start_xd_min(1):-stepRange:0
            fprintf('Current: %d' , xd_min1);
            tic;
            for xd_min2 = start_xd_min(2):-stepRange:0
                temp_res = zeros(start_xd_min(3), 1);
                idx = 1:stepRange:start_xd_min(3);
                parfor (xd_min3 = 1:numel(idx), 16)
                    b = a;
                    temp_xd_min = xd_min;
                    temp_xd_min(1) = xd_min1;
                    temp_xd_min(2) = xd_min2;
                    temp_xd_min(3) = idx(xd_min3);

                    b = simulate(b, temp_xd_min);

                    if (b.isCrashed == true)
                        continue
                    end
                    temp_res(xd_min3) = b.HC;
                end
                if (isempty(temp_res(temp_res>0)))
                    continue
                end
                if (best_HC > min(best_HC, min(temp_res(temp_res>0))))
                    best_HC = min(best_HC, min(temp_res(temp_res>0)));
                end
            end
            best_HC
            fprintf(' in %d seconds\n' , toc);
        end
    case char('rmhca')
        stop = 0;
        for k = 1:generations
            q = round(rand(1)*(n-1)) + 1;
            param = round(rand(1)*stepRange) + 1;
            xd_min(q) = best_xd_min(q) - param;
            a = simulate(a, xd_min);
            if (a.isCrashed ~= true && a.HC < best_HC)
                best_HC = a.HC;
                best_xd_min = a.xd(1:n);
                last_good_xdmin = xd_min;
                stop = 0;
            else
                stop = stop + 1;
               xd_min = last_good_xdmin; 
            end
            if stop == learningBreak
                generationsDone = k;
                break;
            end
        end
    case char('ga')
        individuals = zeros(n,generationSize);
        
        best_set = zeros(n,1);
        best_fitness = 0;
        stop = 0;
        
        for individualNo = 1:generationSize
            for node = 1:n
                individuals(node, individualNo) = round(rand(1)*start_xd_min(node));
            end
        end
        
        gaProcessArchive = GAProcessArchive(topologyFilename)

        while generationsDone < generations
            generationsDone = generationsDone + 1;
            
            gaProcessArchive.bestHCCourse(generationsDone) = best_HC;
            gaProcessArchive.bestFitnessCourse(generationsDone) = best_fitness;
            
            individualFitnesses = zeros(1,generationSize);
            individualUsed = zeros(1,generationSize);
            
            unproductivity = 0;
            
            for individualNo = 1:generationSize
                xd_min = individuals(:, individualNo);
                a = simulate(a, xd_min);
                individualFitnesses(individualNo) = a.fitness;
                
                gaProcessArchive.fitnessCourse(:, end+1) = [generationsDone;a.fitness];
                gaProcessArchive.HCCourse(:, end+1) = [generationsDone;a.HC];
                
                if individualFitnesses(individualNo) > best_fitness
                   best_fitness = individualFitnesses(individualNo);
                   best_xd_min = xd_min
                   best_HC = a.HC
                   resultGenerationNo = generationsDone
                   
                   gaProcessArchive.bestHCFixes(:, end+1) = [generationsDone;best_HC];
                   gaProcessArchive.bestFitnessFixes(:, end+1) = [generationsDone;best_fitness];
                   
                   stop = 0;
                else
                    unproductivity = unproductivity + 1;
                    if unproductivity == generationSize
                       stop = stop + 1; 
                    end
                end
            end
            
            pie((individualFitnesses-0.99) / (sum(individualFitnesses)-generationSize))
            drawnow
            
            if stop == learningBreak
                break;
            end
            
            fitSum = sum(individualFitnesses);
            pairs = [];
            
            while true
                fitRandom = rand;
                
                for individualNo = 1:generationSize
                    if fitRandom < (sum(individualFitnesses(1:individualNo)) / fitSum)
                        if (individualUsed(individualNo) == 1)
                            break;
                        end
                        pairs = [pairs, individuals(:,individualNo)];
                        individualUsed(individualNo) = 1;
                        break;
                    end
                end
                
                if size(pairs, 2) == generationSize
                    break;
                elseif size(pairs, 2) == generationSize-1
                    for individualNo = 1:individualNo
                        if (individualUsed(individualNo) == 0)
                            pairs = [pairs, individuals(:,individualNo)];
                            individualUsed(individualNo) = 1;
                            break;
                        end
                    end
                end
            end
            
            i = round(rand(1)*(n-1)) + 1;
            
            for index = 1:2:generationSize
                individuals(1:i,index) = pairs(1:i,index);
                individuals(i+1:end,index) = pairs(i+1:end,index+1);
                individuals(1:i,index+1) = pairs(1:i,index+1);
                individuals(i+1:end,index+1) = pairs(i+1:end,index);
            end
            
            for individualNo = 1:generationSize
                for node = 1:n
                    randomMutation = rand;
                    if randomMutation < mutationProbability
                        individuals(node, individualNo) = round(rand(1)*start_xd_min(node));
                    end
                end
            end
        end
    otherwise
        disp('Optimization method not found! Try again!\n')
end

a = simulate(a, best_xd_min);

TimeSpent = toc(beginApplicationTic);

% DEBUG HOLDING COST
a.xd(1:n)
round(initial_HC)
round(a.HC)
a.satisfiedRate

TimeSpent;

if fileGenerator
    directory = 'reports';
    filename = sprintf('report_%d_%d_', simTime, m);
    extension = 'txt';
    
    filepath = [directory '/' filename methodType datestr(now, '_mmdd_HHMMSS') '.' extension];
    fid = fopen(filepath,'w');
    
    print_initial_HC = round(initial_HC);
    print_best_HC = round(best_HC);
    print_initial_mf = ((inv((eye(n) + temp)*inv(Lambda)) * (start_xd_min-1)) - dmean) ./ dstd;
    print_best_mf = ((inv((eye(n) + temp)*inv(Lambda)) * (best_xd_min-1)) - dmean) ./ dstd;
    print_satisfied_rate = a.satisfiedRate;
    
    if ~isempty(topologyFilename)
        fprintf(fid, [topologyFilename '\n\n']);
    end
    
    fprintf(fid, 'Generations\n%d/%d (break: %d)\n', generationsDone, generations, learningBreak);
    fprintf(fid, 'Step range\n%d\n', stepRange);
    fprintf(fid, 'Computing time\n%f\n', TimeSpent);
    fprintf(fid, 'Holding cost\n%d -> %d\n', print_initial_HC, print_best_HC);
    fprintf(fid, 'Satisfied rate\n%f\n', print_satisfied_rate);
    if strcmp(methodType, 'ga')
        fprintf(fid, '-> Genetic algorithm parameters\n');
        fprintf(fid, 'Alpha\n%f\n', a.alpha);
        fprintf(fid, 'Beta\n%f\n', a.beta);
        fprintf(fid, 'Fitness\n%f\n', a.fitness);
        fprintf(fid, 'Generation size\n%f\n', generationSize);
        fprintf(fid, 'Mutation probability\n%f\n', mutationProbability);
        fprintf(fid, 'Result generation number\n%f\n', resultGenerationNo);
        
        archiveFilepath = [directory '/' filename methodType datestr(now, '_mmdd_HHMMSS') '.mat'];
        save(archiveFilepath, 'gaProcessArchive');
    end
    fprintf(fid, 'Reference stock levels\nFrom To\n');
    fprintf(fid, '%d\t\t%d \n', [start_xd_min best_xd_min]');
    fprintf(fid, 'Magic factors\nFrom To\n');
    fprintf(fid, '%d\t\t%d \n', [print_initial_mf print_best_mf]');
    fclose(fid);
end

if chartVisibility
    a = simulate(a, [914;523;293]);
    plotter = ChartGenerator(time, simTime, n, a.x, a.u_hist, a.d, a.h);
    stock_level(plotter, 5);
    order_quantity(plotter, 6);
    demand(plotter, 7);
    satisfied_demand(plotter, 8);

    a = simulate(a, start_xd_min);

    plotter = ChartGenerator(time, simTime, n, a.x, a.u_hist, a.d, a.h);
    stock_level(plotter, 1);
    order_quantity(plotter, 2);
    demand(plotter, 3);
    satisfied_demand(plotter, 4);
end
