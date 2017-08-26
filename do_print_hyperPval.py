from RTNet import RTNet
import sys
import os

net = RTNet()
corr_thresh = 0.75

GERM_LAYERS = ['allderm', 'ectoderm', 'mesoderm','endoderm','ESCs-NPCs','ESCs--mesothe', 'ESCs--liver']
outdir = 'pval'
if not os.path.exists(outdir):
	os.makedirs(outdir)

for each_layer in GERM_LAYERS:
	'''
	# Print Hypergeometric P-value of overlap between Neph Switching RT and Neph Switching TRN
	'''
	fout = open(outdir+'/pval_{}.txt'.format(each_layer), 'w')
	net.PrintHyperPVal(which_layer=each_layer, corr_thresh=corr_thresh, fout=fout)
	fout.close()