function pswarm
res=[];reps=1;
fun = @(x)(LightSheetOptimisation(x));
nvars = 20;
if nvars==1
    figure(1);fplot(fun,[0 10])
end
ub = 1*ones(nvars,1);
lb = 0*ub;

% for ex=1:10
sszarr=[5 10 50 100 500 1000 5000];
for ex=6:6
    ssz=sszarr(ex);
    tic
%     ssz=2^ex;


options = optimoptions('particleswarm','SwarmSize',ssz,'Display','off','HybridFcn',@patternsearch);
% rng default  % For reproducibility

soln=zeros(reps,nvars);
for ii=1:reps
[soln(ii,:),fval(ii),exitflag(ii),output] = particleswarm(fun,nvars,lb,ub,options);
itr(ii)=output.iterations;
end
t=toc;
res=[res;ssz/10^3 t mean(soln) std(soln) mean(exitflag) std(exitflag) mean(itr) std(itr)]
save('res.mat','res')
end
res