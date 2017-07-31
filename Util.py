
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
