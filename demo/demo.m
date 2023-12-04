% This is a demo of the model
addpath(genpath('../peripheral'))

b = 1;% distance decay of attractive force
sl = 1;% increment step length L_s

% generate network with default settings. To change number of brain regions and other settings, 
% see documentation of function connectomee_from_pathfinding()
[c_und,node_coord,directed_axons,undirected_axons] = connectomee_from_pathfinding(b,sl,'rng_state',3); % seed for reproduction
n_nodes = size(node_coord,1);

%% evaluate generated network
n_dim = 2; %2-dimensional circle (for drawing)

% weight-distance associations
fs = 50;%fontsize for figure
d = squareform(pdist(node_coord));
indx = find(c_und~=0);
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

%show median axon trajectories
w_max = 1e3;
show_median_axons_specify_max(c_und,undirected_axons,triu(c_und~=0,1),w_max,n_dim);
axis square;axis off;title('');
title('');

% weight distribution (normalized)
c_corrected = c_und ./ sum(c_und,2);%normalize by node strength
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
xlim([-4.5,0]);ylim([0,0.5]);
xticks([-4,0]);yticks([0,0.5]);
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

ks_dim = [.5,.2,.2,.1];
annotation('textbox',ks_dim,'String',sprintf('ks = %.2f',ksstat),'FitBoxToText','on','FontSize',50,'FontWeight','bold','EdgeColor','none');
