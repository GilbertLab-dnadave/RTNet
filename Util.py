
import operator as op

class Util(object):
	def __init__(self):
		pass
	@staticmethod
	def ncr(n, r):
		r = min(r, n-r)
		if r == 0: return 1
		numer = reduce(op.mul, xrange(n, n-r, -1), 1)
		denom = reduce(op.mul, xrange(1, r+1), 1)
		return numer//denom
	@staticmethod
	def GetP(K, i, M, m):
		return Util.ncr(K, i) * Util.ncr(M-K, m-i)
	@staticmethod
	def GetColumn(wb, sheet_name, whichCol, num_max_row):
		ret_list = []
		sheet = wb[sheet_name]
		for row in xrange(2, num_max_row+1):
			val = sheet[whichCol+str(row)]
			ret_list.append(val.value.replace(u'\u2212', u'-'))
		return ret_list