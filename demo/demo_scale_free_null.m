% This is a demo of evaluating the scale free property of the null model
% 300 nodes, threshold at 5% density 
% The following code may take hours to run

c1 = fix(clock);
addpath(genpath('../peripheral_naive')) 

rng(3)%for reproducibility

% same set of parameters as in the model
% beta not needed because growth directions are sampled randomly
pert = 1;
sl = 1;

n_dim = 100;
r = 30;
max_step = 3 .* r ./ sl;
n_nodes = 300;
n_axons = 2e5;
delta_th = pi/12;

return_directed = false;
return_undirected = false;

% generate null networks and perform scale-free test
n_test = 2;
eval_den = 0.05;
all_p = zeros(n_test,1);
all_kmin = zeros(n_test,1);
net_density = zeros(n_test,1);
th_val = zeros(n_test,1);
all_k = zeros(n_test,n_nodes);
all_alpha = zeros(n_test,1);
for i = 1:n_test
    node_rad = sample_circle_rad(n_nodes,pert);%evenly sample nodes followed by perturbation
    node_coord = rad2xy(node_rad,r);%convert polar coordinate to cartesian

    axon_rad = sample_circle_rad(n_axons,2); %random sample axons, 2nd argument > 1 for random sampling
    axon_coord = rad2xy(axon_rad,r);%convert polar coordinate to cartesian

    %simulate the axons: axons - trajectory of each axon; success - if the axon succeed
    [axons,success] = simulate_random_walk_network(r,axon_coord,sl,max_step,delta_th);

    %convert axons to connectivity matrix:
    %c - directed connectivity matrix;
    [c,~,~] = axons2c(axons,success,node_coord,return_directed,return_undirected);

    %makae undirected connectivity matrix
    c_und = c + c' - 2 * diag(diag(c));
    net_density(i) = nnz(c_und)/n_nodes/(n_nodes-1);

    % threshold and test for scale free
    c_th = threshold_proportional(c_und,eval_den);
    th_val(i) = min(c_th(c_th~=0));
    k = sum(c_th~=0);
    all_k(i,:) = k;
    [alpha, kmin, L] = plfit(k,'range',[1.001:0.001:12.001],'finite');
    [p,gof]=plpva(k, kmin,'range',alpha,'reps',1000,'silent');
    all_p(i) = p;
    all_alpha(i) = alpha;
    all_kmin(i) = kmin;
    disp(i);
    disp(median(all_p(1:i)));
end
c2 = fix(clock);
save('scale_free_null.mat','all_p','all_k','all_alpha','all_kmin','net_density','th_val');

%% plot media p degree distribution

[~,select] = min(abs(all_p - median(all_p)));
select_p = all_p(select);
x = all_k(select,:)';
th = min(eval_den,sum(x)/n_nodes/(n_nodes-1));

alpha = all_alpha(select);
xmin = all_kmin(select);

n = length(x);
q = unique(x);
c = hist(x,q)'./n;

c = [[q; q(end)+1] 1-[0; cumsum(c)]]; c(c(:,2)<10^-10,:) = [];
cf = ((xmin:q(end))'.^-alpha)./(zeta(alpha) - sum((1:xmin-1).^-alpha));
cf = [(xmin:q(end)+1)' 1-[0; cumsum(cf)]];
cf(:,2) = cf(:,2) .* c(c(:,1)==xmin,2);

edges = 1:(size(all_k,2)-1);

% ER random binomial/Poisson

% binomial
x_val = edges;
n_trial = n_nodes - 1;
p1 = th;
px = 1- binocdf(x_val,n_trial,p1);
figure;plot(x_val,px,'b','linewidth',5);

xlim([10^0,10^2]);ylim([10^-3,10^0]);
xticks([10^0,10^2]);yticks([10^-3,10^0]);

hold;plot(c(:,1),c(:,2),'r','LineWidth',5);
plot(cf(:,1),cf(:,2),'linewidth',10,'Color','k','LineStyle',':');
set(gca,'Xscale','log','Yscale','log');
Annotate(gca,'textbox',[10^1.2,10^2],[10^-0.4,10^0],'backgroundcolor', 'none', 'string', sprintf('p = %.2f',select_p),'fontsize',40,'fontweight','bold','edgecolor','none');
% Annotate(gca,'textbox',[1.3,2],[-0.4,0],'backgroundcolor', 'none', 'string', sprintf('p = %.2f',select_p),'fontsize',35,'fontweight','bold','edgecolor','none');
legend({'ER','model',sprintf('y ~ k^{%.2f}',-alpha)},'fontsize',40,'location','southwest');
set(gca,'fontsize',40,'fontweight','bold')