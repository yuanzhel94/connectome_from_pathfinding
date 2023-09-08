function nets = random_nets_bu(n_nodes,n_edges,n_nets)
%function nets = random_nets_bu(n_nodes,n_edges)
%   This function generate a number of random networks (n_nets) with
%   specified number of nodes (n_nodes) and undirected edges (n_edges).

nets = zeros(n_nodes,n_nodes,n_nets);
triu_indx = find(triu(ones(n_nodes),1));
n_edge_avail = length(triu_indx);

if n_edge_avail < n_edges
    error("network with %d nodes can have %d undirected edges at maximum, %d edges are asked to generate");
end

for i = 1:n_nets
    indx = randsample(triu_indx,n_edges);
    neti = zeros(n_nodes);
    neti(indx) = 1;
    neti = neti + neti';
    nets(:,:,i) = neti;
end

end