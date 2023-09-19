# connectome_from_pathfinding

## Intro

## Guideline for use
###Key functions
- Peripheral_naive folder contains functions for implementing the model; key functions are listed below. See usage examples in demos.
- sample_circle_rad: sample node center / axons coordinates on the circle, in radian
- rad2xy: convert sampled radian to cartesian
- simulate_network: simulate axons with specified guiding rules
- simulate_random_walk_network: simulated axons with specified random walk rules
- axons2c: assign axons to nodes to form connectivity matrix (directed), need to be converted to undirected

###Demos
- demo.m: generate a model network and test for negative associations between weights and distance, and lognormal weight distribution.
- demo_null.m: generate a random walk null model and test for negative associations between weights and distance, and lognormal weight distribution.
- demo_scale_free.m: evaluate the scale free of generated networks.
- demo_scale_free_null.m: evaluate the scale free of null networks.
- demo_topology.m: evaluate the topology (CC, CPL, SW, and modularity Q) of a generated network. 
- Due to high computational cost, we did not provide a demo for parameter optimization; however, guidance for parameter optimization can be found in demo_topology.m

## Abstract
