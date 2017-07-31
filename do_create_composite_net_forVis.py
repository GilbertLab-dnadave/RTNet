from RTNet import RTNet

net = RTNet()
which_layer = "mesoderm"
corr_thresh = 0.98


'''
# Create Composite Net (RT + reduced TRN) for visualization
# it will return None, but creates csv file for the input for the Cytoscape.
'''
net.CreateCompositeNetForVis(which_layer=which_layer, corr_thresh=corr_thresh)
