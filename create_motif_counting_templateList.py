import networkx as nx
import itertools
import networkx.algorithms.isomorphism as iso
import pydot
import pickle
import sys
import os

def main(argv):
	GList=[]
	n=int(argv[0])
	outdir_pickle = 'DATA_motif_finding_template_list'
	if not os.path.exists(outdir_pickle):
		os.makedirs(outdir_pickle)
	if os.path.exists(outdir_pickle + "/{}nodes_template.pickle".format(n)):
		sys.exit(1)

	numPermu = len(list(itertools.permutations(range(n),2)))
	numCombi = len(list(itertools.combinations(range(n),2)))
	numLoops = n
	# numCasesTotal = numPermu + numCombi + n
	numCasesTotal = numPermu + numCombi

	# it should be different based on 'n'
	if n==2:
		edgeList = [ (0,1), (1,0), (0,1)]
	elif n==3:
		edgeList = [ (0,1), (1,0), (1,2), (2,1), (0,2), (2,0), (0,1), (1,2), (0,2) ]
	elif n==4:
		edgeList= [ (0,1), (1,0), (0,2), (2,0), (0,3), (3,0), (1,2), (2,1), (1,3), (3,1), (2,3), (3,2), (0,1), (0,2), (0,3), (1,2), (1,3), (2,3) ]
	else:
		sys.exit(1)

	for numSeledge in range(1,numCasesTotal+1):
		totalCombiSeledge = itertools.combinations(range(numCasesTotal),numSeledge)
		for eachCase in totalCombiSeledge:
			hereG = nx.MultiDiGraph()
			hereOnlyRTG = nx.Graph()
			for smallerEachCase in eachCase:
				if smallerEachCase<numPermu:
					# default color is red
					hereG.add_edge(*edgeList[smallerEachCase])
				elif numPermu<= smallerEachCase <numPermu + numCombi:
					nodeIndices = edgeList[smallerEachCase]
					swappedNodeIndices = (nodeIndices[1], nodeIndices[0])
					hereG.add_edge(*nodeIndices,color='black')
					hereG.add_edge(*swappedNodeIndices,color='black')
					hereOnlyRTG.add_edge(*nodeIndices)
				else:
					# hereG.add_edge(*edgeList[smallerEachCase])
					pass
			if len(hereG.nodes())!=n:
				totalNodes=range(n)
				for nowNode in hereG.nodes():
					totalNodes.remove(nowNode)
				hereG.add_nodes_from(totalNodes)
			if len(hereOnlyRTG.nodes())!=n:
				totalNodes=range(n)
				for nowNode in hereOnlyRTG.nodes():
					totalNodes.remove(nowNode)
				hereOnlyRTG.add_nodes_from(totalNodes)
			if nx.is_connected(hereOnlyRTG):
				GList.append(hereG)

	duplicatedCase = set()
	em = iso.categorical_multiedge_match('color','red')
	for rawI in range(len(GList)):
		if rawI not in duplicatedCase:
			for rawJ in range(rawI+1,len(GList)):
				if rawJ not in duplicatedCase:
					if nx.is_isomorphic(GList[rawI], GList[rawJ], edge_match=em):
						duplicatedCase.add(rawJ)
	notDuplicatedCase = set(range(len(GList)))
	notDuplicatedCase -= duplicatedCase
	notDuplicatedCase = list(notDuplicatedCase)
	notDuplicatedCase.sort()

	selectedGList = [i for rawI, i in enumerate(GList) if rawI in notDuplicatedCase]

	
	pickle.dump(selectedGList, open(outdir_pickle + '/{}nodes_template.pickle'.format(n),'w'))

	outdir = 'MotifCountingTemplates/{}nodes_template'.format(n)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	templates = pickle.load(open(outdir_pickle + '/{}nodes_template.pickle'.format(n),'r'))
	for rawI, i in enumerate(templates):
		# extracted undirected graph
		extractedUndir = []
		for eachEdge in i.edges(data=True):
			if 'color' in eachEdge[2]:
				extractedUndir.append(eachEdge)
		extractedUndirGraph = nx.Graph()
		extractedUndirGraph.add_edges_from(extractedUndir)
		extractedDir = []
		for eachEdge in i.edges(data=True):
			if 'color' not in eachEdge[2]:
				extractedDir.append((eachEdge[0], eachEdge[1]))
		###
		hereG = nx.MultiDiGraph()
		for j in extractedUndirGraph.edges():
			hereG.add_edge(j[0], j[1], color='black', arrowhead='none')
		for j in extractedDir:
			hereG.add_edge(*j, color='red')
		graphOutName=outdir + '/{}nodes_{}.pdf'.format(n, rawI)
		dotName=outdir + '/{}nodes_{}.dot'.format(n, rawI)
		nx.nx_agraph.write_dot(hereG,dotName)
		graph = pydot.graph_from_dot_file(dotName)
		if type(graph) == type([]):
			graph[0].write_pdf(graphOutName, prog='neato')
		else:
			graph.write_pdf(graphOutName, prog='neato')

if __name__ == '__main__':
	main(sys.argv[1:])
