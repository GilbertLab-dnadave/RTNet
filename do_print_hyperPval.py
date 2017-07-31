from RTNet import RTNet

net = RTNet()
which_layer = "mesoderm"
corr_thresh = 0.98

'''
# Print Hypergeometric P-value of overlap between Neph Switching RT and Neph Switching TRN
'''
net.PrintHyperPVal(which_layer=which_layer, corr_thresh=corr_thresh)
