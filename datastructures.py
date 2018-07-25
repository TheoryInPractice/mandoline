class Bimap:
	def __init__(self):
		self.to = {}
		self.fro = {}

	def __contains__(self, u):
		return u in self.to or u in self.fro

	def __getitem__(self, u):
		if u in self.to:
			return self.to[u]
		return self.fro[u]

	def items(self):
		return self.to.items()

	def put(self, u, v):
		if u in self:
			raise KeyError("{} already in bimap".format(u))
		if v in self:
			raise KeyError("{} already in bimap".format(v))
		self.to[u] = v
		self.fro[v] = u

	def put_all(self, pairs):
		for u,v in pairs:
			self.put(u,v)

	def __repr__(self):
		return ','.join(map(str, self.items()))