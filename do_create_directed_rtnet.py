from RTNet import RTNet

net = RTNet()
which_layer = "ESCs--liver"
etolOpt = 'ltoe'
source = 'APOB'

'''
# Construction of Directed RT network from one source node
# it will return None, but creates csv file for the input for the Cytoscape.
# uncomment to run this version
'''
net.CreateDirectedRTNetFromOneSource(which_layer=which_layer, etolOpt=etolOpt, source=source)

'''
# Construction of Directed RT network from all possible source nodes
# it will return None, but creates csv file for the input for the Cytoscape.
# uncomment to run this version
'''
#net.CreateDirectedRTNetFromAllSources(which_layer=which_layer, etolOpt=etolOpt)

'''
# Construction of Directed RT network from one source recursively
# For example, construction of Drected RTNet from one source to other, and call recursively from other to next other.
# it will return None, but creates csv file for the input for the Cytoscape.
# uncomment to run this version
'''
# net.CreateDirectedRTNetFromOneSourceRecur(which_layer=which_layer, etolOpt=etolOpt, source=source)
