# connectome_from_pathfinding
- This repository contains code for the study: "A generative model of the connectome with dynamic axon growth". 
- A novel connectome generative model is proposed to generate macroscopic connectomes from microscopic, chemoaffinity guided axon outgrowth.

>### Abstract
>Coming soon


>### Guideline to use the model
>- Function ***peripheral/connectome_from_pathfinding.m*** could be used to generate a connectome with specified force decay parameter ***$\beta$*** and step length parameter ***$L_s$***. 
>- $\beta$ and $L_s$ are the two mandatory inputs to the function. ***Unless specified, default implementation of the model will be used for other model parameters (e.g., number of nodes, number of axons, etc). These parameters can be specified with "name-value" pairs.*** See function peripheral/connectome_from_pathfinding.m for more details.

>### Subdirectories
>- #### peripheral
>   This subdir contains the functions required to replicate major findings in the paper, including:
>   -  connectome_from_pathfinding.m - one-line code for generating connectoems from specified parameters.
>   -  supporting functions to connectome_from_pathfinding.m
>   -  functions that generate random walk null model
>   -  functions that visualize generated axons
>   -  third-party codes for analysis and figure plots, including Brain connectivity toolbox, power-law fit, and annotation.
>   -  and more features ...
>- #### demo
>   This subdir contains a series of demo that replicates major findings in the paper, including:
>   -  ***demo.m*** - replicates axon visualization, weight-distance associations, and weight distributions for the ***model networks***.
>   -  ***demo_null.m*** - replicate the same sets of properties as in demo.m, but for random walk null model networks.
>   -  ***demo_scalefree.m*** - replicate the scale-free degree distribution for model networks.
>   -  ***demo_scalefree_null.m*** - replicate the scale-free evaluation for random walk null networks.
>   -  ***demo_topology.m*** - demonstrate how the normalized complex topological measures (i.e., CC, CPL, SW, Q) were evaluated and optimized.


