
from RTNet import RTNet

net = RTNet()
which_layer = "mesoderm"
corr_thresh = 0.98

'''
# Construction of Bipartite Network between
# gene expression network with all 3 fold difference genes and RT network with all switching genes
# it will return None, but creates csv file for the input for the Cytoscape.
'''
net.CreateBinet(which_layer, corr_thresh)
