# Thu Nov 19 14:22:52 2020
# Author: Jeffrey Durieux, MSc

# What: simulationscript for simulation 1 of clusterwise joint ica


#### Design ####

# 20 repetitions of each combination of the following factors:
#   - number of Q, 2, 5, 10
#   - number of R, 2, 3, 4
#   - N per cluster 20, 30, 50
#   - one-to-one correlation: no, low (.25), medium (.50)
#   - SNR 19, 4, .33 (5%, 40%, 75%)
# 
# Fixed: 5000 voxels, 100 random starts

#### correlated clusters
