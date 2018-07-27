from helpers import short_str

class Bimap:
    """
        A bijection between two sets. Anything not explicilty
        put into this index structures is assumed to be paired
        with itself.
    """
    def __init__(self):
        self.to = {}
        self.fro = {}

    def __contains__(self, u):
        return u in self.to or u in self.fro

    def __getitem__(self, u):
        if u in self.to:
            return self.to[u]
        if u in self.fro:
            return self.fro[u]
        return u

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
        return ','.join(map(short_str, self.items()))


class Indexmap:
    """
        Maps a set of vertices (e.g. arbitrary, hashable objects)
        to indices [0,...,n-1].
    """
    def __init__(self, size):
        self.vertex_to_index = {}
        self.index_to_vertex = [None] * size

    def __len__(self):
        return len(self.index_to_vertex)

    def put(self, index, vertex):
        if index >= len(self) or index < 0:
            raise IndexError()
        self.vertex_to_index[vertex] = index
        self.index_to_vertex[index] = vertex

    def vertex_at(self, index):
        if index >= len(self) or index < 0 or self.index_to_vertex[index] == None:
            raise IndexError()
        return self.index_to_vertex[index]

    def index_of(self, vertex):
        return self.vertex_to_index[vertex]

    def indices_of(self, vertices):
        return [self.vertex_to_index[v] if v in self.vertex_to_index else None for v in vertices]

    def vertices_at(self, indices):
        return [self.index_to_vertex[i] for i in indices]

    def __repr__(self):
        return 'IM['+','.join(map(str,self.index_to_vertex))+']'        