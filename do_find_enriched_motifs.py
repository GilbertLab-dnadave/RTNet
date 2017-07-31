from RTNet import RTNet

net = RTNet()
which_layer = "mesoderm"
corr_thresh = 0.98


'''
# Find Enriched motifs in Composite Net (RT + reduced TRN)
# it will return None, but creates pdf file.
'''
net.FindEnrichedMotifs(which_layer=which_layer, corr_thresh=corr_thresh)
