
'''
# Construction of Bipartite Network (with specific selected EXP genes) between
# gene expression network that consists of selected EXP genes and RT network with all switching genes
# it will return None, but creates csv file for the input for the Cytoscape.
'''

from RTNet import RTNet
r = RTNet()

GERM_LAYERS = ['ESCs--liver', 'ESCs--panc']
CORR_THRESH_LIST = [0.75, 0.9]
SHEETS = ["TopGenes", "Only TFs"]
COLS = ['E', 'H'] # Liver_D16, and Panc_D12
MAX_ROW = {'TopGenes': [100, 101], 'Only TFs':[50, 51]}
RATIOS = [0.5, 0.75]
for each_col_ind, (each_layer, each_col) in enumerate(zip(GERM_LAYERS, COLS)):
	for each_corr in CORR_THRESH_LIST:
		for each_sheet in SHEETS:
			for each_ratio in RATIOS:
				r.CreateBinetWithSpecificSelectedEXPGenes(each_layer, each_corr, each_sheet, each_col, MAX_ROW[each_sheet][each_col_ind], each_ratio)
