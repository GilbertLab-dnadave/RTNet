###### Required version of Python
Python 2.7

###### Install dependencies:
```
	> pip install -r requirements.txt
```

##### Create community RT network For RT Networks visualization. It will return None, but creates txt file for the input for the SAFE algorithm.
```
	> python do_create_community_rtnet_forVis.py
```

##### Construction of Directed RT network. It will return None, but creates csv file for the input for the Cytoscape.
```
	> python do_create_directed_rtnet.py
```

##### Print Hypergeometric P-value of overlap between Neph Switching RT and Neph Switching TRN
```
	> python do_print_hyperPval.py
```

##### Create Composite Net (RT + reduced TRN) for visualization. It will return None, but creates csv file for the input for the Cytoscape.
```
	> python do_create_composite_net_forVis.py
```

##### Find Enriched motifs in Composite Net (RT + reduced TRN). It will return None, but creates pdf file.
```
	> python do_find_enriched_motifs.py
```

##### Construction of Bipartite Network between gene expression network with all 3 fold difference genes and RT network with all switching genes. It will return None, but creates csv file for the input for the Cytoscape.
```
	> python do_create_binet.py
```
