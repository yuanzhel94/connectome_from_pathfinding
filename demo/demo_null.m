% This is a demo of the random walk null model
addpath(genpath('../peripheral'))

rng(3);%seed for reproduction

% same set of parameters as in the model
% growth direction is sampled randomly rather than computed from forces
pert = 1;
sl = 1;

n_dim = 2;
r = 30;
max_step = 3 .* r ./ sl;
n_nodes = 84;
n_axons = 2e5;
self_attract = true;
theta = pi/12;

return_directed = true;
return_undirected = true;

% sample nodes and axons: same as the model
node_rad = sample_circle_rad(n_nodes,pert);
node_coord = rad2xy(node_rad,r);
axon_rad = sample_circle_rad(n_axons,2); % random sample axons
axon_coord = rad2xy(axon_rad,r);


% simulate random walk null
[axons,success] = simulate_random_walk_network(r,axon_coord,sl,max_step,theta);
[c,directed_axons,undirected_axons] = axons2c(axons,success,node_coord,return_directed,return_undirected);
c_und = c + c' - 2 * diag(diag(c));

%% evaluate generated null network
fs = 50; % font size for plots
d = squareform(pdist(node_coord))
indx = find(c_und~=0);

% weight distance association
figure;scatter(d(indx),log10(c_und(indx)),30,[1,1,1] .* 0.5,'filled');
h = lsline;
set(h,'linewidth',5,'color','k');
set(gca,'fontsize',fs,'fontweight','bold');
rval = corr(d(indx),log10(c_und(indx)));
text_dim = [.4,.8,.2,.1];
annotation('textbox',text_dim,'String',sprintf('r = %.2f',rval),'FitBoxToText','on','FontSize',fs,'FontWeight','bold','EdgeColor','none');
xticks([0,60]);
xlim([0,60]);ylim([0,4]);
xticks([0,60]);yticks([0,4]);



% show all generated connections (despite too tense for visualization, see below)
w_max = 1e3;
show_median_axons_specify_max(c_und,undirected_axons,triu(c_und~=0,1),w_max,n_dim);
axis square;axis off;
title('');

% show 1000 randomly selected axons
show_axons(axons(1:1000,:),'k',1);
axis square;axis off;title('');

% null network is too dense for visualization. Visualize 5% connections only
triu_indx = find(triu(ones(size(c_und)),1));
mask_indx = randsample(triu_indx,round(length(triu_indx)*0.05));
mask_mat = zeros(size(c_und));
mask_mat(mask_indx) = 1;
mask_mat = mask_mat + mask_mat';
c_mask = c_und .* mask_mat;
show_median_axons_specify_max(c_mask,undirected_axons,triu(c_mask~=0,1),w_max,n_dim);
axis square;axis off;
title('');


% weight distribution
c_corrected = c_und ./ sum(c_und,2);
off_diag = find(ones(n_nodes) - diag(diag(ones(n_nodes))));
hist_vals = log10(c_corrected(off_diag));
hist_vals = hist_vals(~isinf(hist_vals));
[vals,bin] = histcounts(hist_vals,15);

fig_inset = figure;histogram(hist_vals,bin,'normalization','pdf','EdgeAlpha',0,'FaceColor','r');
set(gca,'fontsize',fs,'FontWeight','bold');
fitted = fitdist(hist_vals,'Normal');
hold;
y = pdf(fitted,bin);
plot(bin,y,'k','linewidth',5)
[h,p,ksstat,cv] = kstest(hist_vals,fitted);
%     ks_dim = [.25,.2,.2,.1];
%     annotation('textbox',ks_dim,'String',sprintf('ks stat = %.2f',ksstat),'FitBoxToText','on','FontSize',fs,'FontWeight','bold','EdgeColor','none');
xlim([-4.5,0]);ylim([0,4]);
xticks([-4,0]);yticks([0,4]);

xval = sort(unique(hist_vals),'ascend');
cdf_norm = cdf(fitted,xval);
fig_main = figure;hold;
rep = 1000;
R_cdf = zeros(rep,length(xval));
for j = 1:rep
    Ri = random(fitted,size(hist_vals));
    [cdf_Ri,xj] = ecdf(Ri);
    x_diff = xval' - xj;
    x_eval = xval' >= xj;
    cdf_eval = repmat(cdf_Ri,1,length(xval));
    cdf_eval(~x_eval) = 0;
    R_cdf(j,:) = max(cdf_eval,[],1);
end
R_5prc = prctile(R_cdf,5,1);
R_95prc = prctile(R_cdf,95,1);
patch([xval;flipud(xval)],[R_5prc,fliplr(R_95prc)],'k');
mdl_cdf = cdfplot(hist_vals);
%     plot(xval,cdf_norm,'LineWidth',5,'color','k');
xlabel('');ylabel('');set(gca,'fontsize',fs,'fontweight','bold');
set(mdl_cdf,'LineWidth',5,'color','r');

grid off;
box off;
title('');
xlim([-4.5,0]);ylim([0,1]);
xticks([-4,0]);yticks([0,1]);

[hm,hi] = inset(fig_main,fig_inset);
set(hi,'position',[0.25,0.6,0.25,0.25 ],'fontsize',30)

ks_dim = [.54,.2,.2,.1];
annotation('textbox',ks_dim,'String',sprintf('ks = %.2f',ksstat),'FitBoxToText','on','FontSize',50,'FontWeight','bold','EdgeColor','none');
