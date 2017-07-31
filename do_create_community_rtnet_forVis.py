from RTNet import RTNet

net = RTNet()
which_layer = "mesoderm"
corr_thresh = 0.98
degree_filter_thresh = 20

'''
# Create community RT network For RT Networks visualization
# it will return None, but creates txt file for the input for the SAFE algorithm.
'''
net.CreateCommunityRTNetForVis(which_layer=which_layer, corr_thresh=corr_thresh, degree_filter_thresh=degree_filter_thresh)
