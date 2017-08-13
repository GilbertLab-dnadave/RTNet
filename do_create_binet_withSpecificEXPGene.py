from RTNet import RTNet
r = RTNet()

GERM_LAYERS = ['ESCs--mesothe', 'ESCs--smooth', 'ESCs--liver', 'ESCs--panc']
SPECIFIC_EXP_GENE = ['FOXA1', 'PDX1']
for each_layer in GERM_LAYERS:
	for each_gene in SPECIFIC_EXP_GENE:
		r.CreateBinetWithSpecificEXPGene(each_layer, 0.9, each_gene, 0.75)
