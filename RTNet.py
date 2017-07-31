'''
Dependency:
	numpy : pip install numpy
	networkx : pip install networkx
	python-louvain (community) : pip install python-louvain
	pydot : pip install pydot
'''

import itertools
from operator import itemgetter
import sys
import networkx as nx
import os
import numpy as np
import community
from collections import defaultdict
from Util import Util
from math import log, exp
import networkx.algorithms.isomorphism as iso
import pickle
import copy
import random
import scipy.stats as st
from shutil import copyfile
from heapq import heappush
from heapq import heappop

class RTNet(object):
	def __init__(self, source="./uniq_rt_exp.csv", neph_source="./DATA_NephNet/allderm_network.txt"):
		self.GERM_LAYERS = ['ectoderm', 'mesoderm', 'endoderm', 'ESCs-NPCs', 'ESCs-neural-MSCs', 'ESCs--mesothe', 'ESCs--smooth',\
		 'ESCs--liver', 'ESCs--panc']
		self.GERM_LAYERS_MATCHING_IND = [[0, 1, 15, 16, 17, 18], [0, 1, 9, 10, 11, 12, 13, 14], [0, 1, 2, 3,4,5, 6,7,8], \
		[0, 1, 18], [0, 1, 15, 16, 17], \
		[0, 1, 9, 10, 11, 12, 13], [0, 1, 9, 10, 11, 14], \
		[0, 1, 2, 3,4,5], [0, 1, 2, 6,7,8]]
		# read source data
		data = open(source, 'r').readlines()
		data=[i.split(',') for i in data]
		data=[[i.rstrip() for i_ind, i in enumerate(each) if i_ind!=65] for each in data]
		# remove heading row
		data = data[1:]
		# remove a row which has at least one empty data
		complete_data_withsex= []
		for each_line in data:
			if any(i=='' for i in each_line):
				pass
			else:
				complete_data_withsex.append(each_line)

		self.GENE_NAMES_WITHSEX = [itemgetter(1)(each) for each in complete_data_withsex]
		self.EXP_DATA_WITHSEX = [[float(i) for i in itemgetter(*range(8, 27))(each)] for each in complete_data_withsex]
		self.RT_DATA_WITHSEX = [[float(i) for i in itemgetter(*range(27, 46))(each)] for each in complete_data_withsex]
		PRE_LOC_WITHSEX = [list(itemgetter(*([2,5,6,7]))(each)) for each in complete_data_withsex]
		self.LOC_WITHSEX = [list(itemgetter(*[0,2])(each)) if each[1]=='+' else list(itemgetter(*[0, 3])(each)) for each in PRE_LOC_WITHSEX ]
		self.LOC_WITHSEX = [[i if i_ind==0 else int(i) for i_ind, i in enumerate(j)] for j in self.LOC_WITHSEX]
		# remove sex chromosomes
		complete_data = [[col for col_ind, col in enumerate(each)] for each in complete_data_withsex if each[2]!='chrX' and each[2]!='chrY']
		self.GENE_NAMES = [itemgetter(1)(each) for each in complete_data]
		self.EXP_DATA = [[float(i) for i in itemgetter(*range(8, 27))(each)] for each in complete_data]
		self.RT_DATA = [[float(i) for i in itemgetter(*range(27, 46))(each)] for each in complete_data]
		PRE_LOC = [list(itemgetter(*([2,5,6,7]))(each)) for each in complete_data]
		self.LOC = [list(itemgetter(*[0,2])(each)) if each[1]=='+' else list(itemgetter(*[0, 3])(each)) for each in PRE_LOC ]
		self.LOC = [[i if i_ind==0 else int(i) for i_ind, i in enumerate(j)] for j in self.LOC]
		# read source of neph data
		neph = open(neph_source, 'r').readlines()
		neph = [i.split() for i in neph]
		neph = [[i.strip() for i in j] for j in neph]
		neph_gene_set = set()
		for each_edge in neph:
			neph_gene_set.add(each_edge[0])
			neph_gene_set.add(each_edge[1])
		all_gene_set = set(self.GENE_NAMES)
		neph_gene_set = neph_gene_set & all_gene_set
		neph_gene_set = list(neph_gene_set)
		neph_gene_set.sort()
		self.NEPH_GENE_NAMES = neph_gene_set

		# self.NEPH_GENE_NAMES = open('./tmp_geneNames.csv', 'r').readlines()
		# self.NEPH_GENE_NAMES = [i.strip() for i in self.NEPH_GENE_NAMES]

	# which_layer should be string. For example, allderm or ectoderm or ESCs--liver or ....
	# data_type should be string. For example, RT or EXP.
	# return dictionary {"GENE_NAMES": ... , "LOC": ..., "DATA": ...}
	def GetWhichLayerData(self, which_layer, data_type):
		ret_dic = {}
		ret_dic['GENE_NAMES'] = self.GENE_NAMES
		ret_dic['LOC'] = self.LOC
		if which_layer == 'allderm':
			if data_type == 'RT':
				ret_dic['DATA'] = self.RT_DATA
				return ret_dic
			elif data_type == 'EXP':
				ret_dic['DATA'] = self.EXP_DATA
				return ret_dic
			else:
				print "Please input exact data_type string."
				sys.exit(1)
		try:
			that_ind = self.GERM_LAYERS.index(which_layer)
		except ValueError:
			print "Please input exact which_layer string."
			sys.exit(1)
		if data_type == 'RT':
			ret_dic['DATA'] = [list(itemgetter(*self.GERM_LAYERS_MATCHING_IND[that_ind])(each)) for each in self.RT_DATA]
			return ret_dic
		elif data_type == 'EXP':
			ret_dic['DATA'] = [list(itemgetter(*self.GERM_LAYERS_MATCHING_IND[that_ind])(each)) for each in self.EXP_DATA]
			return ret_dic
		else:
			print "Please input exact data_type string."
			sys.exit(1)
	# which_layer should be string. For example, allderm or ectoderm or ESCs--liver or ....
	# data_type should be string. For example, RT or EXP.
	# return dictionary {"GENE_NAMES": ... , "LOC": ..., "DATA": ...}
	def GetWhichLayerDataWithsex(self, which_layer, data_type):
		ret_dic = {}
		ret_dic['GENE_NAMES'] = self.GENE_NAMES_WITHSEX
		ret_dic['LOC'] = self.LOC_WITHSEX
		if which_layer == 'allderm':
			if data_type == 'RT':
				ret_dic['DATA'] = self.RT_DATA_WITHSEX
				return ret_dic
			elif data_type == 'EXP':
				ret_dic['DATA'] = self.EXP_DATA_WITHSEX
				return ret_dic
			else:
				print "Please input exact data_type string."
				sys.exit(1)
		try:
			that_ind = self.GERM_LAYERS.index(which_layer)
		except ValueError:
			print "Please input exact which_layer string."
			sys.exit(1)
		if data_type == 'RT':
			ret_dic['DATA'] = [list(itemgetter(*self.GERM_LAYERS_MATCHING_IND[that_ind])(each)) for each in self.RT_DATA_WITHSEX]
			return ret_dic
		elif data_type == 'EXP':
			ret_dic['DATA'] = [list(itemgetter(*self.GERM_LAYERS_MATCHING_IND[that_ind])(each)) for each in self.EXP_DATA_WITHSEX]
			return ret_dic
		else:
			print "Please input exact data_type string."
			sys.exit(1)
	# check if this gene is switching (RT).
	@staticmethod
	def IsSwitching(rt):
		row = [float(i) for i in rt]
		e_check=0
		l_check=0
		for i in xrange(len(row)):
			if row[i]>=0.3:
				e_check += 1
			if row[i]<=-0.3:
				l_check += 1
		if e_check!=0 and l_check!=0:
			return True
		else:
			return False
	# which_layer should be string. For example, allderm or ectoderm or ESCs--liver or ....
	# return dictionary {"GENE_NAMES": ... , "LOC": ..., "DATA": ...}
	# DATA is always RT type.
	def GetWhichLayerSwitchingData(self, which_layer):
		pre_dic = self.GetWhichLayerData(which_layer, "RT")
		GENE_NAMES = pre_dic['GENE_NAMES']
		LOC = pre_dic['LOC']
		DATA = pre_dic['DATA']
		sel_ind = []
		for each_row_ind, each_row in enumerate(DATA):
			if RTNet.IsSwitching(each_row):
				sel_ind.append(each_row_ind)
		ret_dic = {}
		ret_dic['GENE_NAMES'] = list(itemgetter(*sel_ind)(GENE_NAMES))
		ret_dic['LOC'] = list(itemgetter(*sel_ind)(LOC))
		ret_dic['DATA'] = list(itemgetter(*sel_ind)(DATA))
		return ret_dic
	# which_layer should be string. For example, allderm or ectoderm or ESCs--liver or ....
	# return dictionary {"GENE_NAMES": ... , "LOC": ..., "DATA": ...}
	# DATA is always RT type.
	def GetWhichLayerSwitchingDataWithsex(self, which_layer):
		pre_dic = self.GetWhichLayerDataWithsex(which_layer, "RT")
		GENE_NAMES = pre_dic['GENE_NAMES']
		LOC = pre_dic['LOC']
		DATA = pre_dic['DATA']
		sel_ind = []
		for each_row_ind, each_row in enumerate(DATA):
			if RTNet.IsSwitching(each_row):
				sel_ind.append(each_row_ind)
		ret_dic = {}
		ret_dic['GENE_NAMES'] = list(itemgetter(*sel_ind)(GENE_NAMES))
		ret_dic['LOC'] = list(itemgetter(*sel_ind)(LOC))
		ret_dic['DATA'] = list(itemgetter(*sel_ind)(DATA))
		return ret_dic
	# Construction of RT networks based on coordinated changes in RT
	# it will return networkx graph
	def CreateRTNet(self, which_layer, corr_thresh, bool_only_neph=False, bool_out_txt=False):
		if bool_only_neph is False:
			outdir = 'RTNet/rawRTNet'
			net_dic = self.GetWhichLayerSwitchingData(which_layer)
		else:
			outdir = 'RTNet/rawRTNet_onlyNeph'
			net_dic = self.GetWhichLayerNephSwitchingData(which_layer)

		if bool_out_txt:
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			fout = open(outdir + '/{}_{}_rawRTNet.csv'.format(which_layer, corr_thresh), 'w')
			fout.write("GENE1,GENE2,CORR,DISTANT,DISTANCE\n")

		GENE_NAMES = net_dic['GENE_NAMES']
		LOC = net_dic['LOC']
		DATA = net_dic['DATA']
		G = nx.Graph()

		for i in GENE_NAMES:
			G.add_node(i)

		for ind_i in xrange(len(DATA)):
			for ind_j in xrange(ind_i+1, len(DATA)):
				here_corr = np.corrcoef(DATA[ind_i], DATA[ind_j])[0,1]
				if here_corr >= corr_thresh:
					distance = -1
					if LOC[ind_i][0] == LOC[ind_j][0]:
						distance = abs(LOC[ind_i][1] - LOC[ind_j][1])
					distant = 1
					if LOC[ind_i][0] == LOC[ind_j][0] and distance < 500000:
						distant = 0
					if distant==1:
						'''
						commonPartitionNum, partitionCount, aPartitionNum, and bPartitionNum are
							for community network for the visualization of RTNet
						commonPartitionNum : if gene1 and gene2 for an edge in RTNet are in same partition,
							commonPartitionNum will have some value rather than -1.
						partitionCount : total number of edges that belongs to this commonPartitionNum
						aPartitionNum : partition number for gene1
						bPartitionNum : partition number for gene2
						'''
						G.add_edge(GENE_NAMES[ind_i], GENE_NAMES[ind_j], corr = here_corr, \
						commonPartitionNum = -1, partitionCount = 0, aPartitionNum = -1, bPartitionNum = -1)
						if bool_out_txt:
							out_str = "{},{},{},{},{}\n".format(GENE_NAMES[ind_i], GENE_NAMES[ind_j], here_corr, distant, distance)
							fout.write(out_str)
		return G
	# Create community RT network For RT Networks visualization
	# it will return None, but creates txt file for the input for the SAFE algorithm.
	def CreateCommunityRTNetForVis(self, which_layer, corr_thresh, degree_filter_thresh):
		outdir = 'RTNet/visRTNet/{}_{}'.format(which_layer, corr_thresh)
		G = self.CreateRTNet(which_layer, corr_thresh)

		if not os.path.exists(outdir):
			os.makedirs(outdir)

		fout_com = open(outdir + '/{}_{}_{}_community_RTNet.csv'.format(which_layer, corr_thresh, degree_filter_thresh), 'w')
		fout_com.write("GENE1,GENE2,CORR,COMMON_PARTITION_NUM,COMMON_PARTITION_COUNT,A_PARTITION_NUM,B_PARTITION_NUM\n")
		fout_num_count = open(outdir + '/{}_{}_partitionNumCountInfo.csv'.format(which_layer, corr_thresh), 'w')
		fout_num_count.write("COMMON_PARTITION_NUM,COMMON_PARTITION_COUNT\n")
		fout_gsafe = open(outdir + '/{}_{}_{}_graphForSAFE.txt'.format(which_layer, corr_thresh, degree_filter_thresh), 'w')

		partition = community.best_partition(G)
		partitionNumList = list(set(partition.values()))
		partitionNumList.sort()

		# assign commonPartitionNum, particionCount, aPartitionNum, and bPartitionNum
		partitionNumCount = defaultdict(int)
		for eachEdge in G.edges():
			eachEdgeDict = G[eachEdge[0]][eachEdge[1]]
			eachEdgeDict['aPartitionNum'] = partition.get(eachEdge[0])
			eachEdgeDict['bPartitionNum'] = partition.get(eachEdge[1])
			G.node[eachEdge[0]]['partitionNum'] = partition.get(eachEdge[0])
			G.node[eachEdge[1]]['partitionNum'] = partition.get(eachEdge[1])
			if partition.get(eachEdge[0]) == partition.get(eachEdge[1]):
				eachEdgeDict['commonPartitionNum'] = partition.get(eachEdge[0])
			parNum = eachEdgeDict['commonPartitionNum']
			if parNum!=-1:
				partitionNumCount[parNum] += 1
		for eachEdge in G.edges():
			eachEdgeDict = G[eachEdge[0]][eachEdge[1]]
			if eachEdgeDict['commonPartitionNum'] != -1:
				eachEdgeDict['partitionCount'] = partitionNumCount[eachEdgeDict['commonPartitionNum']]
		#####
		# assign degree of counterpart partition number for all nodes
		for eachNode in G.nodes():
			for eachPartitionNum in partitionNumList:
				nodeDict = G.node[eachNode]
				nodeDict[eachPartitionNum] = 0
			for eachPartitionNum in partitionNumList:
				nodeDict = G.node[eachNode]
				for eachTargetNode in G[eachNode]:
					if G.node[eachTargetNode]['partitionNum'] == eachPartitionNum:
						nodeDict[eachPartitionNum] += 1
		#####

		# after sorting, write file for community RT Network.
		sortedEdges = sorted(G.edges(data=True), key=lambda x:(x[2]['partitionCount'], x[2]['commonPartitionNum']), reverse=True)
		for eachEdge in sortedEdges:
			eachEdgeDict = G[eachEdge[0]][eachEdge[1]]
			degree1 = G.node[eachEdge[0]][eachEdgeDict['bPartitionNum']]
			degree2 = G.node[eachEdge[1]][eachEdgeDict['aPartitionNum']]
			if degree1 > degree_filter_thresh and degree2 > degree_filter_thresh:
				out_str = '{},{},{},{},{},{},{}\n'.format(eachEdge[0], eachEdge[1], eachEdgeDict['corr'],\
					eachEdgeDict['commonPartitionNum'], eachEdgeDict['partitionCount'],\
					eachEdgeDict['aPartitionNum'], eachEdgeDict['bPartitionNum'])
				fout_com.write(out_str)
				fout_gsafe.write("{}\t{}\t{}\t{}\t{}\n".format(eachEdge[0], eachEdge[0], eachEdge[1], eachEdge[1], 1))
		#####

		# write file for common partition num count info.
		for partitionNum, partitionCount in sorted(partitionNumCount.items(), key=lambda x:x[1], reverse=True):
			out_str = '{},{}\n'.format(partitionNum, partitionCount)
			fout_num_count.write(out_str)
	# which_layer should be string. For example, ESCs--liver, ESC--mesothe or ....
	# return dictionary {"GENE_NAMES": ... , "LOC": ..., "DATA": ...}
	# DATA is always RT type.
	def GetWhichLayerSwitchingDataForDirected(self, which_layer):
		GERM_LAYERS = ["ESCs-neural-MSCs", "ESCs--mesothe", "ESCs--smooth","ESCs--liver","ESCs--panc"]
		dupIndexSet=[ [1, 3], [1,4,6], [1,4], [1], [1]]
		try:
			that_ind = GERM_LAYERS.index(which_layer)
		except ValueError:
			print "Please input exact which_layer string."
			sys.exit(1)

		pre_dic = self.GetWhichLayerData(which_layer, "RT")
		GENE_NAMES = pre_dic['GENE_NAMES']
		LOC = pre_dic['LOC']
		DATA = pre_dic['DATA']

		dupIndex = dupIndexSet[that_ind]
		selIndex = list(set(range(len(DATA[0]))) - set(dupIndex))
		selIndex.sort()

		# ESC second removed and dup removed
		DATA = [list(itemgetter(*selIndex)(each)) for each in DATA]

		sel_ind = []
		for each_row_ind, each_row in enumerate(DATA):
			if RTNet.IsSwitching(each_row):
				sel_ind.append(each_row_ind)
		ret_dic = {}
		ret_dic['GENE_NAMES'] = list(itemgetter(*sel_ind)(GENE_NAMES))
		ret_dic['LOC'] = list(itemgetter(*sel_ind)(LOC))
		ret_dic['DATA'] = list(itemgetter(*sel_ind)(DATA))
		return ret_dic
	'''
	it checks if rtEL is straight or not.
	for example, straight is L-L-E-E or E-E-E-L
	'''
	@staticmethod
	def HelperStraight(firstNum, mNum, thirdNum, etolOpt, totalLen, rtEL):
		if etolOpt == 'ltoe':
			firstCheck = 'L'
			thirdCheck = 'E'
		elif etolOpt == 'etol':
			firstCheck = 'E'
			thirdCheck = 'L'
		template = [firstCheck] * firstNum + ['M'] * mNum + [thirdCheck] * thirdNum
		if template == rtEL:
			return True
		return False
	'''
	it checks if rtEL is a sequence of particular first and second and whatever after.
	for example, if firstNum is 1 and secondNums is 1 and etolOpt is ltoe, then L-E-...
	'''
	@staticmethod
	def HelperWhatever(firstNum, secondNum, etolOpt, totalLen, rtEL):
		if etolOpt == 'ltoe':
			firstCheck = 'L'
			secondCheck = 'E'
		elif etolOpt == 'etol':
			firstCheck = 'E'
			secondCheck = 'L'
		additional = [secondCheck] * (totalLen-firstNum-secondNum)
		template = [firstCheck] * firstNum + [secondCheck] * secondNum + additional
		if template==rtEL: return True
		for i in xrange(1, totalLen-firstNum-secondNum+1):
			additional[-i] = firstCheck
			template = [firstCheck] * firstNum + [secondCheck] * secondNum + additional
			if template == rtEL: return True
		return False
	def HelperInitDirectedRTNet(self, which_layer, etolOpt):
		if etolOpt not in ['etol', 'ltoe']:
			print "Please input exact etolOpt string."
			sys.exit(1)

		pre_dic = self.GetWhichLayerSwitchingDataForDirected(which_layer)
		GENE_NAMES = pre_dic['GENE_NAMES']
		DATA = pre_dic['DATA']
		rtEL = [['E' if i>=0.3 else 'L' if i<=-0.3 else 'M' for i in each] for each in DATA]
		return {"rtEL": rtEL, 'GENE_NAMES': GENE_NAMES}

	# Construction of Directed RT network from one source node
	# it will return None, but creates csv file for the input for the Cytoscape.
	def CreateDirectedRTNetFromOneSource(self, which_layer, etolOpt, source, bool_out_csv=True, fout_param=None):
		if bool_out_csv is False and fout_param is None or bool_out_csv is True and fout_param is not None:
			print "Please use legal input combinations!"
			sys.exit(1)

		pre_dic = self.HelperInitDirectedRTNet(which_layer, etolOpt)
		rtEL = pre_dic["rtEL"]
		GENE_NAMES = pre_dic["GENE_NAMES"]

		try:
			source_ind = GENE_NAMES.index(source)
		except ValueError:
			print "Source {} is not possible.".format(source)
			return

		if bool_out_csv:
			outdir = 'RTNet/DirectedRTNet/FromOneSource'
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			fout = open(outdir + "/{}_{}_{}.csv".format(which_layer, etolOpt, source), 'w')
			fout.write("GENE1,GENE2,GENE1_Hierarchy,GENE2_Hierarchy\n")
		else:
			fout = fout_param

		for target_ind in xrange(len(rtEL)):
			if source_ind == target_ind: continue
			fEL = rtEL[source_ind]; sEL = rtEL[target_ind]
			for tl2 in xrange(1, len(fEL)-1):
				if RTNet.HelperWhatever(tl2, 1, etolOpt, len(fEL), fEL):
					# one step diff
					if RTNet.HelperStraight(tl2+1, 0, len(fEL)-tl2-1, etolOpt, len(sEL), sEL):
						# hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(fEL) +','+str(sEL) + '\n'
						hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(len(fEL)-tl2) +','+str(len(fEL)-tl2-1) + '\n'
						fout.write(hereStr)
					# two step diff
					elif tl2!=len(fEL)-2 and RTNet.HelperStraight(tl2+2, 0, len(fEL)-tl2-2, etolOpt, len(sEL), sEL):
						# hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(fEL) +','+str(sEL) + '\n'
						hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(len(fEL)-tl2) +','+str(len(fEL)-tl2-2) + '\n'
						fout.write(hereStr)
				# m case
				if tl2!=len(fEL)-2 and RTNet.HelperWhatever(tl2, 2, etolOpt, len(fEL), fEL) and RTNet.HelperStraight(tl2+1, 1, len(sEL)-tl2-2, etolOpt, len(sEL), sEL):
					# hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(fEL) +','+str(sEL) + '\n'
					hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(len(fEL)-tl2) +','+str(len(fEL)-tl2-1) + '\n'
					fout.write(hereStr)
	# Construction of Directed RT network from all possible source nodes
	# it will return None, but creates csv file for the input for the Cytoscape.
	def CreateDirectedRTNetFromAllSources(self, which_layer, etolOpt):
		pre_dic = self.HelperInitDirectedRTNet(which_layer, etolOpt)
		rtEL = pre_dic["rtEL"]
		GENE_NAMES = pre_dic["GENE_NAMES"]

		outdir = 'RTNet/DirectedRTNet/FromAllSources'
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		fout = open(outdir + "/{}_{}.csv".format(which_layer, etolOpt), 'w')
		fout.write("GENE1,GENE2,GENE1_Hierarchy,GENE2_Hierarchy\n")
		for source_ind, each_gene in enumerate(GENE_NAMES):
			self.CreateDirectedRTNetFromOneSource(which_layer, etolOpt, each_gene, False, fout)
	# Construction of Directed RT network from one source recursively
	# For example, construction of Drected RTNet from one source to other, and call recursively from other to next other.
	# it will return None, but creates csv file for the input for the Cytoscape.
	def CreateDirectedRTNetFromOneSourceRecur(self, which_layer, etolOpt, source):
		pre_dic = self.HelperInitDirectedRTNet(which_layer, etolOpt)
		rtEL = pre_dic["rtEL"]
		GENE_NAMES = pre_dic["GENE_NAMES"]

		try:
			source_ind = GENE_NAMES.index(source)
		except ValueError:
			print "Source {} is not possible.".format(source)
			return

		outdir = 'RTNet/DirectedRTNet/FromOneSourceRecur'
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		fout = open(outdir + "/{}_{}_{}.csv".format(which_layer, etolOpt, source), 'w')
		fout.write("GENE1,GENE2,GENE1_Hierarchy,GENE2_Hierarchy\n")

		self.HelperCreateDirectedRTNetFromOneSourceRecur(which_layer, etolOpt, source, fout)

	def HelperCreateDirectedRTNetFromOneSourceRecur(self, which_layer, etolOpt, source, fout):
		pre_dic = self.HelperInitDirectedRTNet(which_layer, etolOpt)
		rtEL = pre_dic["rtEL"]
		GENE_NAMES = pre_dic["GENE_NAMES"]

		try:
			source_ind = GENE_NAMES.index(source)
		except ValueError:
			print "Source {} is not possible.".format(source)
			return

		for target_ind in xrange(len(rtEL)):
			if source_ind == target_ind: continue
			fEL = rtEL[source_ind]; sEL = rtEL[target_ind]
			for tl2 in xrange(1, len(fEL)-1):
				if RTNet.HelperWhatever(tl2, 1, etolOpt, len(fEL), fEL):
					# one step diff
					if RTNet.HelperStraight(tl2+1, 0, len(fEL)-tl2-1, etolOpt, len(sEL), sEL):
						# hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(fEL) +','+str(sEL) + '\n'
						hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(len(fEL)-tl2) +','+str(len(fEL)-tl2-1) + '\n'
						fout.write(hereStr)
						self.HelperCreateDirectedRTNetFromOneSourceRecur(which_layer, etolOpt, GENE_NAMES[target_ind], fout)
					# two step diff
					elif tl2!=len(fEL)-2 and RTNet.HelperStraight(tl2+2, 0, len(fEL)-tl2-2, etolOpt, len(sEL), sEL):
						# hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(fEL) +','+str(sEL) + '\n'
						hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(len(fEL)-tl2) +','+str(len(fEL)-tl2-2) + '\n'
						fout.write(hereStr)
						self.HelperCreateDirectedRTNetFromOneSourceRecur(which_layer, etolOpt, GENE_NAMES[target_ind], fout)
				# m case
				if tl2!=len(fEL)-2 and RTNet.HelperWhatever(tl2, 2, etolOpt, len(fEL), fEL) and RTNet.HelperStraight(tl2+1, 1, len(sEL)-tl2-2, etolOpt, len(sEL), sEL):
					# hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(fEL) +','+str(sEL) + '\n'
					hereStr = GENE_NAMES[source_ind] + ',' + GENE_NAMES[target_ind]+ ','+ str(len(fEL)-tl2) +','+str(len(fEL)-tl2-1) + '\n'
					fout.write(hereStr)
					self.HelperCreateDirectedRTNetFromOneSourceRecur(which_layer, etolOpt, GENE_NAMES[target_ind], fout)
	# which_layer should be string. For example, allderm or ectoderm or ESCs--liver, ESC--mesothe or ....
	# return dictionary {"GENE_NAMES": ... , "LOC": ..., "DATA": ...}
	def GetWhichLayerNephSwitchingData(self, which_layer):
		pre_dic = self.GetWhichLayerData(which_layer, "RT")
		GENE_NAMES = pre_dic['GENE_NAMES']
		LOC = pre_dic['LOC']
		DATA = pre_dic['DATA']
		# only choosing neph genes
		sel_ind = []
		for each_name in self.NEPH_GENE_NAMES:
			here_ind = GENE_NAMES.index(each_name)
			sel_ind.append(here_ind)
		GENE_NAMES = list(itemgetter(*sel_ind)(GENE_NAMES))
		LOC = list(itemgetter(*sel_ind)(LOC))
		DATA = list(itemgetter(*sel_ind)(DATA))
		#####
		# filter out non-switching genes
		sel_ind = []
		for each_row_ind, each_row in enumerate(DATA):
			if RTNet.IsSwitching(each_row):
				sel_ind.append(each_row_ind)
		ret_dic = {}
		ret_dic['GENE_NAMES'] = list(itemgetter(*sel_ind)(GENE_NAMES))
		ret_dic['LOC'] = list(itemgetter(*sel_ind)(LOC))
		ret_dic['DATA'] = list(itemgetter(*sel_ind)(DATA))
		return ret_dic
	# it will return networkx DiGraph
	def ReadNephTRNNet(self, which_layer):
		GERM_LAYERS = ['allderm', 'ectoderm', 'mesoderm', 'endoderm', 'ESCs-NPCs', 'ESCs--mesothe', 'ESCs--liver']
		if which_layer not in GERM_LAYERS:
			print "Please input exact which_layer string."
			sys.exit(1)
		NephG=nx.DiGraph()
		NephNet = open('./DATA_NephNet/{}_network.txt'.format(which_layer)).readlines()
		NephNet = [i.split() for i in NephNet]
		for each_edge in NephNet:
			NephG.add_edge(each_edge[0], each_edge[1])
		return NephG
	# Print Hypergeometric P-value of overlap between Neph Switching RT and Neph Switching TRN
	def PrintHyperPVal(self, which_layer, corr_thresh):
		NephG = self.ReadNephTRNNet(which_layer)

		rt_dic = self.GetWhichLayerNephSwitchingData(which_layer)
		RT_NAMES = rt_dic["GENE_NAMES"]

		print "Number of nodes and edges of the Neph graph"
		print len(NephG.nodes())
		print len(NephG.edges())
		print

		n = len(NephG.nodes())
		m = len(NephG.edges())
		M = 2 * Util.ncr(n, 2)

		# only neph RT
		RTG = self.CreateRTNet(which_layer, corr_thresh, True)
		# try composing RT network (RT network is undirected.)

		NephG_edgeDict={key: 0 for key,value in zip(NephG.edges(), xrange(len(NephG.edges())))}
		NephG_edgeSet = set(NephG_edgeDict.keys())
		RTG_edgeSet = set([(a, b) for a, b in RTG.edges()])

		print "Number of nodes and edges of the RTG graph"
		print len(RTG.nodes())
		print len(RTG.edges())
		print

		K = len(RTG.edges())

		common_edges = 0
		for each_edge_rt in RTG.edges():
			if each_edge_rt in NephG_edgeDict or (each_edge_rt[1], each_edge_rt[0]) in NephG_edgeDict:
				common_edges += 1
		rt_edges = 0
		neph_edges = 0
		neph_edges = len(NephG_edgeSet) - common_edges
		rt_edges = len(RTG_edgeSet) - common_edges
		print "Number of common edges between Neph and RT when direction is not significant"
		print common_edges
		print
		k2 = common_edges
		print "Only Neph Edges and Only RT Edges"
		print neph_edges
		print rt_edges
		print
		# Get PVal of hypergeo (direction is not significant).
		res = 0
		for i in xrange(k2+1, K+1):
			res += Util.GetP(K, i, M, m)
		res = log(res) - log(Util.ncr(M, m))
		res = exp(res)
		print "When direction is not significant, hypergeometric p-value is"
		print res
		print
	# Get Neph TRN (non-RT nodes are removed) and Neph RT Graph
	def Get_RTNodesNephG_RTG(self, which_layer, corr_thresh):
		NephG = self.ReadNephTRNNet(which_layer)
		# only neph RT
		RTG = self.CreateRTNet(which_layer, corr_thresh, True)

		# reduce Neph G
		willBeRemoved=[]
		rt_node_set = set(RTG.nodes())
		for n in NephG.nodes():
			if n not in rt_node_set:
				willBeRemoved.append(n)
		for w in willBeRemoved:
			NephG.remove_node(w)
		#####
		ret_dic = {}
		ret_dic['NephG'] = NephG
		ret_dic['RTG'] = RTG
		return ret_dic
	# Create Composite Net (RT + reduced TRN) for visualization
	# it will return None, but creates csv file for the input for the Cytoscape.
	def CreateCompositeNetForVis(self, which_layer, corr_thresh):
		pre_dic = self.Get_RTNodesNephG_RTG(which_layer, corr_thresh)
		NephG = pre_dic['NephG']
		RTG = pre_dic['RTG']

		NephG_edgeDict={key: 0 for (key,value) in zip(NephG.edges(), xrange(len(NephG.edges())))}

		composedGraph=nx.MultiDiGraph()
		for eachEdge in RTG.edges():
			composedGraph.add_edge(eachEdge[0], eachEdge[1], color='black', arrowhead='none')
			# composedGraph.add_edge(eachEdge[1], eachEdge[0], color='black', arrowhead='none')

		reducedNephG = nx.DiGraph()
		num_nodes_list = [2,3,4]

		for numNodes in num_nodes_list:
			for nodeOfSubG in itertools.combinations(RTG.nodes(), numNodes):
				subG = RTG.subgraph(nodeOfSubG)
				if nx.is_connected(subG):
					permuNodes = itertools.permutations(nodeOfSubG, 2)
					for eachEdge in permuNodes:
						if eachEdge in NephG_edgeDict and eachEdge not in reducedNephG.edges():
							reducedNephG.add_edge(*eachEdge)

		red_count = 0
		for eachEdge in reducedNephG.edges():
			composedGraph.add_edge(*eachEdge, color='red')
			red_count+=1

		dashed_count = 0
		permuNodes = itertools.permutations(RTG.nodes(), 2)
		for eachEdge in permuNodes:
			if eachEdge in NephG_edgeDict and eachEdge not in composedGraph.edges():
				composedGraph.add_edge(eachEdge[0], eachEdge[1], color="dashed")
				dashed_count+=1

		print "Number of nodes and edges of the RT graph"
		print len(RTG.nodes())
		print len(RTG.edges())
		print "Number of nodes and edges of the composed graph"
		print len(composedGraph.nodes())
		print len(composedGraph.edges())
		print "Number of NephTRN directed edges within RTN motif and NephTRN directed edges outside RTN motif"
		print red_count, dashed_count

		outdir = 'CompositeNet/visCompositeNet'.format(which_layer, corr_thresh)
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		fout = open(outdir + '/{}_{}_CompositeNet.csv'.format(which_layer, corr_thresh), 'w')
		fout.write("GENE1,GENE2,COLOR\n")

		for eachEdge in composedGraph.edges(data=True):
			here_str = eachEdge[0] + ',' + eachEdge[1] + ','
			if eachEdge[2].get("color") == 'red':
				here_str += "directed"
			elif eachEdge[2].get("color") == 'black':
				here_str += "undirected"
			elif eachEdge[2].get("color") == 'dashed':
				here_str += "dashed"
			here_str += "\n"
			fout.write(here_str)
		fout.close()
	@staticmethod
	def shuffleGraph(network):
		networkCopy=network.copy()
		shuffledNetwork = network.copy()
		shuffledNetwork.clear()
		edgeSize = len(networkCopy.edges())
		while len(networkCopy.edges())>1:
			firstEdge = random.choice(networkCopy.edges())
			networkCopy.remove_edge(*firstEdge)
			infFallCheck=0
			while True:
				secondEdge = random.choice(networkCopy.edges())
				edgeDupCheckSet=set()
				edgeDupCheckSet |= set(firstEdge)
				edgeDupCheckSet |= set(secondEdge)
				if len(edgeDupCheckSet)==4:
					networkCopy.remove_edge(*secondEdge)
					break
				infFallCheck+=1
				if infFallCheck>=100:
					return False
			shuffledNetwork.add_edge(firstEdge[0], secondEdge[1])
			shuffledNetwork.add_edge(secondEdge[0], firstEdge[1])
			if len(shuffledNetwork.edges())%100 == 0:
				# print('{} done'.format(len(shuffledNetwork.edges())))
				pass
		if len(networkCopy.edges())==1:
			shuffledNetwork.add_edge(*networkCopy.edges()[0])
			networkCopy.clear()
		return shuffledNetwork
	# Get Composite Net (RT + reduced TRN) graph as networkx graph
	def GetCompositeGraphForMotifCounting(self, numNodes, NephG, RTG):
		composedGraph=nx.MultiDiGraph()
		for eachEdge in RTG.edges():
			composedGraph.add_edge(eachEdge[0], eachEdge[1], color='black')
			composedGraph.add_edge(eachEdge[1], eachEdge[0], color='black')
		reducedNephG = nx.DiGraph()
		NephG_edgeDict={key: 0 for (key,value) in zip(NephG.edges(), xrange(len(NephG.edges())))}
		for nodeOfSubG in itertools.combinations(RTG.nodes(), numNodes):
			subG = RTG.subgraph(nodeOfSubG)
			if nx.is_connected(subG):
				permuNodes = itertools.permutations(nodeOfSubG, 2)
				for eachEdge in permuNodes:
					if eachEdge in NephG_edgeDict:
						reducedNephG.add_edge(*eachEdge)
		for eachEdge in reducedNephG.edges():
			composedGraph.add_edge(*eachEdge)
		return composedGraph
	# Count motif
	def CountMotifFromCompositeG(self, numNodes, composedGraph):
		em = iso.categorical_multiedge_match('color','red')
		outdir_pickle = 'DATA_motif_finding_template_list'
		# motif counting
		templates=pickle.load(open(outdir_pickle + '/{}nodes_template.pickle'.format(numNodes),'r'))
		foundMotifList=[]
		foundMotifIndList=[]
		allCounting = [0 for i in xrange(len(templates))]
		nodes = composedGraph.nodes()
		for nodeOfSubG in itertools.combinations(nodes,numNodes):
			subG = composedGraph.subgraph(nodeOfSubG)
			# extracted undirected graph
			extractedUndir = []
			for eachEdge in subG.edges(data=True):
				if 'color' in eachEdge[2]:
					extractedUndir.append(eachEdge)
			extractedUndirGraph = nx.Graph()
			extractedUndirGraph.add_nodes_from(subG.nodes())
			extractedUndirGraph.add_edges_from(extractedUndir)
			###
			if not nx.is_connected(extractedUndirGraph):
				continue
			copiedFoundMotifList = copy.copy(foundMotifList)
			found=False
			for rawInd, eachMotif in enumerate(copiedFoundMotifList):
				if nx.faster_could_be_isomorphic(subG, eachMotif):
					if nx.is_isomorphic(subG, eachMotif, edge_match=em):
						allCounting[foundMotifIndList[rawInd]] += 1
						found=True
						break
			if not found:
				for rawI, templ in enumerate(templates):
					if nx.faster_could_be_isomorphic(subG, templ):
						if nx.is_isomorphic(subG, templ, edge_match=em):
							foundMotifList.append(subG)
							foundMotifIndList.append(rawI)
							allCounting[rawI] += 1
							break
		return allCounting
	# Find Enriched motifs in Composite Net (RT + reduced TRN)
	# it will return None, but creates pdf file.
	def FindEnrichedMotifs(self, which_layer, corr_thresh):
		num_nodes_list = [2,3,4]
		# num_nodes_list = [2,3]
		randomize_iteration = 20
		for numNodes in num_nodes_list:
			os.system('python create_motif_counting_templateList.py {}'.format(numNodes))

		pre_dic = self.Get_RTNodesNephG_RTG(which_layer, corr_thresh)
		NephG_orig = pre_dic['NephG']
		RTG = pre_dic['RTG']

		for numNodes in num_nodes_list:

			NephG = NephG_orig.copy()

			composedGraph = self.GetCompositeGraphForMotifCounting(numNodes, NephG, RTG)
			allCounting = self.CountMotifFromCompositeG(numNodes, composedGraph)
			randomAllCounting_list = []
			print "{}nodes motifs".format(numNodes)
			print allCounting

			for rand_iter in xrange(randomize_iteration):
				# Create Random Graph and
				while True:
					shuffledRT = RTNet.shuffleGraph(RTG)
					if shuffledRT:
						break
					else:
						print "OMG! Shuffling edge fails. Repeat."
				NephG = NephG_orig.copy()
				NephG.remove_edges_from(NephG.selfloop_edges())
				while True:
					shuffledNeph = RTNet.shuffleGraph(NephG)
					if shuffledNeph:
						break
					else:
						print "OMG! Shuffling edge fails. Repeat."
				#####

				hereRandomG = self.GetCompositeGraphForMotifCounting(numNodes, shuffledNeph, shuffledRT)
				randomAllCounting = self.CountMotifFromCompositeG(numNodes, hereRandomG)
				randomAllCounting_list.append(randomAllCounting)

				print "random iter {}: {}".format(rand_iter, randomAllCounting)
			outdir = 'CompositeNet/motifCount/{}nodes/{}_{}'.format(numNodes, which_layer, corr_thresh)
			motif_dir = 'MotifCountingTemplates/{}nodes_template'.format(numNodes)
			if not os.path.exists(outdir):
				os.makedirs(outdir)
			# calculate p-value
			randomAllCounting_list = np.array(randomAllCounting_list)
			randomMean = np.mean(randomAllCounting_list, axis=0)
			randomStdev = np.std(randomAllCounting_list, axis=0)
			pval_list = []
			for rawI, eachSource in enumerate(allCounting):
				if randomStdev[rawI]==0 and randomMean[rawI]==0 and eachSource>0:
					hereZ = 1000.
				elif randomStdev[rawI]==0 and randomMean[rawI]==0 and eachSource==0:
					hereZ = float('nan')
				else:
					hereZ = (eachSource - randomMean[rawI]) / randomStdev[rawI]
				pval = 1.0 - st.norm.cdf(hereZ)
				# if hereZ>1.645:
				if hereZ>=2.54:
					pval_list.append((rawI, pval))
			pval_list.sort(key=lambda x:x[1])
			print pval_list
			for each_ind, each_pval in pval_list:
				copyfile(motif_dir + '/{}nodes_{}.pdf'.format(numNodes, each_ind), \
					outdir + '/pval{:.7f}_motifID{}.pdf'.format(each_pval, each_ind))
			print
			#####
	# Construction of Bipartite Network between
	# gene expression network with all 3 fold difference genes and RT network with all switching genes
	# it will return None, but creates csv file for the input for the Cytoscape.
	def CreateBinet(self, which_layer, corr_thresh):
		nFoldDiff = 3
		pre_dic1 = self.GetWhichLayerData(which_layer, 'EXP')
		EXP_GENE_NAMES = pre_dic1['GENE_NAMES']
		EXP_LOC = pre_dic1['LOC']
		EXP_DATA = pre_dic1['DATA']

		# at least 3-fold difference
		include_exp_ind = []
		for each_exp_ind, each_exp_data in enumerate(EXP_DATA):
			if max(each_exp_data) / min(each_exp_data) > nFoldDiff:
				include_exp_ind.append(each_exp_ind)
		EXP_GENE_NAMES = list(itemgetter(*include_exp_ind)(EXP_GENE_NAMES))
		EXP_LOC = list(itemgetter(*include_exp_ind)(EXP_LOC))
		EXP_DATA = list(itemgetter(*include_exp_ind)(EXP_DATA))
		#####

		pre_dic2 = self.GetWhichLayerSwitchingData(which_layer)
		RT_GENE_NAMES = pre_dic2['GENE_NAMES']
		# append '-rt' at the end of the name
		RT_GENE_NAMES = [i+'-rt' for i in RT_GENE_NAMES]
		RT_LOC = pre_dic2['LOC']
		RT_DATA = pre_dic2['DATA']

		# # only ltoe constraint in RT
		# ltoe_ind = []
		# for each_rt_ind, each_rt in enumerate(RT_DATA):
		# 	first_encountered = None
		# 	for each in each_rt:
		# 		if first_encountered is None and each<=-0.3:
		# 			first_encountered = 'L'
		# 			break
		# 		elif first_encountered is None and each>=0.3:
		# 			first_encountered = 'E'
		# 			break
		# 	if first_encountered == 'L':
		# 		ltoe_ind.append(each_rt_ind)
		# RT_GENE_NAMES = list(itemgetter(*ltoe_ind)(RT_GENE_NAMES))
		# RT_LOC = list(itemgetter(*ltoe_ind)(RT_LOC))
		# RT_DATA = list(itemgetter(*ltoe_ind)(RT_DATA))
		# #####
		G = nx.Graph()
		for exp_ind, each_exp in enumerate(EXP_DATA):
			for rt_ind, each_rt in enumerate(RT_DATA):
				hereCorr = np.corrcoef(each_exp, each_rt)[0,1]
				if hereCorr >= corr_thresh:
					distance = -1
					if RT_LOC[rt_ind][0] == EXP_LOC[exp_ind][0]:
						distance = abs(RT_LOC[rt_ind][1] - EXP_LOC[exp_ind][1])
					distant = 1
					if RT_LOC[rt_ind][0] == EXP_LOC[exp_ind][0] and distance < 500000:
						distant = 0
					if distant==1:
						G.add_edge(EXP_GENE_NAMES[exp_ind], RT_GENE_NAMES[rt_ind], corr = hereCorr)
		# to accelerate run time of finding bipartite net, make complete net only if they are in a connected component.
		Gcopy = G.copy()
		for com in nx.connected_components(Gcopy):
			subG = G.subgraph(com)
			for node in subG.nodes():
				for node2 in subG.nodes():
					if node!=node2:
						if node.endswith('-rt') and node2.endswith('-rt'):
							G.add_edge(node, node2)
						if not (node.endswith('-rt') or node2.endswith('-rt')):
							G.add_edge(node, node2)
		#####
		# get top 20 largest cliques (bipartite net)
		# only maintain 20 largest expLen * rtLen cliques
		cliques=nx.find_cliques(G)
		sortedByLenEdgesTargetCList = []
		for c in cliques:
			rtLen = len([k for k in c if k.endswith('-rt')])
			expLen = len(c) - rtLen
			if rtLen>3 and expLen>3:
				hereTuple = (expLen * rtLen, c)
				if len(sortedByLenEdgesTargetCList)>=20:
					if min(sortedByLenEdgesTargetCList)[0] < rtLen * expLen:
						heappop(sortedByLenEdgesTargetCList)
						heappush(sortedByLenEdgesTargetCList, hereTuple)
				else:
					heappush(sortedByLenEdgesTargetCList, hereTuple)
		sortedByLenEdgesTargetCList = sorted(sortedByLenEdgesTargetCList, key=lambda x: x[0], reverse=True)
		#####

		outdir = "BiNet/{}_{}".format(which_layer, corr_thresh)
		if not os.path.exists(outdir):
			os.makedirs(outdir)

		fileSuffix = 0
		binet_list = []
		for cliqueInd in range(len(sortedByLenEdgesTargetCList)):
			if fileSuffix>=20:
				break
			selectedRt = [k for k in sortedByLenEdgesTargetCList[cliqueInd][1] if k.endswith('-rt')]
			selectedExp = list( set(sortedByLenEdgesTargetCList[cliqueInd][1]) - set(selectedRt) )
			if len(selectedRt) >3 and len(selectedExp)>3:
				selectedExpScore = []
				selectedRtScore = []
				for eachExp in selectedExp:
					avgCorr = 0
					for eachRt in selectedRt:
						avgCorr += G.edge[eachExp][eachRt]['corr']
					avgCorr /= len(selectedRt)
					selectedExpScore.append(avgCorr)
				for eachRt in selectedRt:
					avgCorr = 0
					for eachExp in selectedExp:
						avgCorr += G.edge[eachRt][eachExp]['corr']
					avgCorr /= len(selectedExp)
					selectedRtScore.append(avgCorr)

				binet_list.append((fileSuffix, selectedExp, selectedExpScore, selectedRt, selectedRtScore))
				fileSuffix+=1
		binet_list.sort(key=lambda x:(sum(x[2]) / float(len(x[2])) + sum(x[4]) / float(len(x[4]))), reverse = True)
		for each_ind, each in enumerate(binet_list):
			avgCorr = ( sum(each[2]) / float(len(each[2])) + sum(each[4]) / float(len(each[4])) ) / 2.0
			outFileN = open(outdir + '/{}_{}_{}_{:.6f}.csv'.format(which_layer, corr_thresh, each_ind, avgCorr), 'w')
			outFileN.write("EXP genes:\n")
			for i_ind, i in enumerate(each[1]):
				outFileN.write(i+","+str(each[2][i_ind]) + '\n')
			outFileN.write("RT genes:\n")
			for i_ind, i in enumerate(each[3]):
				outFileN.write(i+","+str(each[4][i_ind]) + '\n')

if __name__ == '__main__':
	argv = sys.argv[1:]
	x = RTNet()
	# z = x.FindEnrichedMotifs(argv[0], float(argv[1]))
	# z = x.CreateCompositeNetForVis(argv[0], float(argv[1]))
	z = x.CreateBinet(argv[0], float(argv[1]))
