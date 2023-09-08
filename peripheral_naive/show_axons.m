function f = show_axons(axons,c,w)
%function f = show_axons(axons)
%   This function plots the trajectory of all input axons

f = figure;
hold;
n_axons = size(axons,1);
for i = 1:n_axons
    axoni = cat(1,axons{i,:});
    plot(axoni(:,1),axoni(:,2),'Color',c,'LineWidth',w);
end