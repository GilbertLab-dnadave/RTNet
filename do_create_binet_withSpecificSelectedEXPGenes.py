
'''
# Construction of Bipartite Network (with specific selected EXP genes) between
# gene expression network that consists of selected EXP genes and RT network with all switching genes
# it will return None, but creates csv file for the input for the Cytoscape.
'''

from RTNet import RTNet
r = RTNet()

GERM_LAYERS = ['ESCs--liver', 'ESCs--panc', 'ESCs--smooth', 'ESCs-neural-MSCs', 'allderm', 'allderm']
COLS = ['E', 'H', 'L', 'O', 'H', 'A'] # Liver_D16, Panc_D12, SM, MSC, Panc_D12, hESC
MAX_ROW = {'TopGenes': [100, 101, 101, 101, 101, 94], 'Only TFs':[50, 51, 51, 51, 51, 44]}
CORR_THRESH_LIST = [0.75, 0.9]
SHEETS = ["TopGenes", "Only TFs"]
RATIOS = [0.5, 0.75]
for each_col_ind, (each_layer, each_col) in enumerate(zip(GERM_LAYERS, COLS)):
	for each_corr in CORR_THRESH_LIST:
		for each_sheet in SHEETS:
			for each_ratio in RATIOS:
				r.CreateBinetWithSpecificSelectedEXPGenes(each_layer, each_corr, each_sheet, each_col, MAX_ROW[each_sheet][each_col_ind], each_ratio)
