clear all; clc; close all;rng(11,'v5normal')
%Original code:https://uk.mathworks.com/matlabcentral/fileexchange/30660-simple-example-of-pso-algorithm

numInd = 40;
numIter = 500; tolerance = [1e-6 1e-6 1e-6 1e-6 1e-6];
range_min=0; % Range for initial swarm's elements
range_max=1;
numVar=20; % Number of variables
iter=1; % Number of iteration
k=1; % weight of stocastic element

 % Vector of swarm's velocity
c_cost = 0;
gamma = 1;
radius1=1000;
radius = [1000 1000 1000 1000 1000];% Initial radius for stop criterion
%PSO
cfn=1;
while cfn<=5
    ind = range_min + (range_max-range_min).*rand(numInd,20);
    v=zeros(numInd,numVar);
    while iter<numIter && radius(cfn)>tolerance(cfn)
        tic
        parfor l=1:numInd
            a = LightSheetOptimisation((ind(l,:))',0);
            valF(l,1)=a(cfn); % Fitness function for the swarm
        end
        [valF_ord,index]=sort(valF); % Sort the objective function's values for the swarm and identify the leader
        leader=ind(index(1),:);
        cst = LightSheetOptimisation(leader',1);
        cost(iter) = cst(cfn);
        parfor l=1:size(ind,1) % Calculates the new velocity and positions for all swarm's elements
            fi=rand();
            v(l,:)=(1-(sqrt(k*fi))/2)*v(l,:)+k*fi*(leader-ind(l,:)); % Velocity
            ind(l,:)=ind(l,:)+gamma*(1-(sqrt(k*fi))/2)*v(l,:)+(1-k*fi)*(leader-ind(l,:)); % Position
        end
        radius(cfn)=norm(leader-ind(index(end),:)); % Calculates the new radius
        fprintf('Iteration Number:%d\t Radius:%f\t Cost:%f\n',iter,radius(cfn),cost(iter));
        if iter>5 && (cost(iter)-cost(iter-1) < 1e-2)
            c_cost = c_cost+1;
        end
        iter=iter+1; % Increases the number of iteration
        
        toc
    end
    filename = ['costfn' num2str(cfn) '.mat'];
    save(filename)
    cfn = cfn+1;
end

p_min=ind(1:20,:); % Output variables
f_min=valF_ord(1:20,:);