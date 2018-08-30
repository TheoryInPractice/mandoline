#!/usr/bin/env python3
from itertools import combinations, chain

def CheckExt(choices):
    import argparse, os.path

    """
        Argparse 'type' to check for specific file extensions.
        See: https://stackoverflow.com/questions/15203829/python-argparse-file-extension-checking
    """
    class Act(argparse.Action):
        def __call__(self,parser,namespace,fname,option_string=None):
            ext = os.path.splitext(fname)[1][1:]
            print(">>",fname, ext)
            if ext not in choices:
                option_string = '({})'.format(option_string) if option_string else ''
                raise argparse.ArgumentError(self, "file doesn't end with one of {}{}".format(choices,option_string))
            else:
                setattr(namespace,self.dest,fname)
    return Act

def suborder(order, selection):
    """ Returns the suborder induced by 'selection' """
    return tuple(filter(lambda s: s in selection, order))

def singlehash(x):
    return (hash(x) * 14695981039346656037) % (1 << 64)

def pairhash(x, y):
    fnv_prime = 1099511628211
    fnv_offset = 14695981039346656037
    modulo = 1 << 64 # 2**64

    res = fnv_offset
    res ^= hash(x)
    res = (res * fnv_prime) % modulo
    res ^= hash(y)
    res = (res * fnv_prime) % modulo

    return res

def vertex_str(l):
    return short_str(encode_vertices(l))

def short_str(l):
    l = list(l)
    if len(l) == 0:
        return '.'
    return ''.join(map(str,l))

def encode_vertices(l):
    l = list(l)
    return [encode_vertex(x) for x in l]

def decode_vertices(s):
    if s == '.':
        return []
    return [decode_vertex(c) for c in s]

def encode_vertex(x):
    assert x >= 0
    if x < 10:
        return str(x) # 0-9
    if x >= 10 and x <= 35:
        return chr(97 + (x-10)) # a-z
    # If necessary: extend to uppercase.
    assert False, "Cannot encode number {}".format(x)

def decode_vertex(c):
    c = ord(c)
    if c >= 48 and c <= 57:
        return c - 48 # 0-9
    if c >= 97 and c <= 122:
        return c - 87 # a-z
    assert False, "Cannot decode character {}".format(c)

if __name__ == "__main__":
    for i in range(100):
        print(i, inthash(inthash(inthash(i))))


def powerset(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(0, len(s)+1))

def powerset_nonempty(iterable):
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(1, len(s)+1))
