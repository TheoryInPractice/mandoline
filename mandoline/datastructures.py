import unittest
from helpers import short_str, pairhash

class Interval:
    def __init__(self, low, high):
        self.low = low
        self.high = high

    def __len__(self):
        return self.high - self.low + 1

    def __contains__(self, x):
        return self.low <= x and x <= self.high

    def __iter__(self):
        return iter(range(self.low, self.high+1))

    def __repr__(self):
        return "[{};{}]".format(self.low, self.high)

class Bimap:
    """
        A bijection between two sets. Anything not explicilty
        put into this index structures is assumed to be paired
        with itself.
    """
    def __init__(self):
        self.to = {}
        self.fro = {}
        self.hash = 14695981039346656037

    def __contains__(self, u):
        return u in self.to or u in self.fro

    def __len__(self):
        return len(self.to)

    def __getitem__(self, u):
        if u in self.to:
            return self.to[u]
        if u in self.fro:
            return self.fro[u]
        return u

    def __hash__(self):
        return self.hash

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        if self.hash != other.hash:
            return False
        if self.to.items() != other.to.items():
            return False
        return True

    def items(self):
        return self.to.items()

    def source(self):
        return set(self.to.keys())

    def target(self):
        return set(self.fro.keys())

    def put(self, u, v):
        if self[u] == v:
            return # Already paired
        if self[u] != u:
            raise KeyError("{} already in bimap".format(u))
        if self[v] != v:
            raise KeyError("{} already in bimap".format(v))
        if u == v:
            return
        self.to[u] = v
        self.fro[v] = u
        self.hash ^= pairhash(u,v)

    def put_all(self, pairs):
        for u,v in pairs:
            self.put(u,v)

    def is_identity(self):
        return len(self.to) == 0

    def __repr__(self):
        return ','.join(map(short_str, self.items()))


class Indexmap:
    """
        Maps a set of vertices (e.g. arbitrary, hashable objects)
        to indices [0,...,n-1].
        TODO: Find better name since it now does not _only_ store vertices.
    """
    def __init__(self, size):
        self.vertex_to_index = {}
        self.index_to_vertex = [None] * size

    def __len__(self):
        return len(self.index_to_vertex)

    def __iter__(self):
        return iter(range(len(self.index_to_vertex)))

    def __getitem__(self, x):
        return self.index_to_vertex[x]

    def order(self):
        return iter(self.index_to_vertex)

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


class TestBimap(unittest.TestCase):

    def test_hashing_equals(self):
        mpA = Bimap()
        mpB = Bimap()
        self.assertEqual(mpA.hash, mpB.hash)
        self.assertEqual(hash(mpA), hash(mpB))
        self.assertEqual(mpA, mpB)

        mpA.put_all([(0,1),(2,3),(4,5),(6,7),(100,200)])
        self.assertNotEqual(mpA, mpB)

        mpB.put_all(reversed([(0,1),(2,3),(4,5),(6,7)]))
        self.assertNotEqual(mpA, mpB)
        mpB.put(100,200)

        self.assertEqual(mpA, mpB)
        mpB.put(10,20)

        self.assertNotEqual(mpA, mpB)
        mpA.put(20,10)
        self.assertNotEqual(mpA, mpB)


if __name__ == '__main__':
    unittest.main()
