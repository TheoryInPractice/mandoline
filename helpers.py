#!/usr/bin/env python3

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


def short_str(l):
    l = list(l)
    if len(l) == 0:
        return '.'
    return ''.join(map(str,l))


if __name__ == "__main__":
    for i in range(100):
        print(i, inthash(inthash(inthash(i))))
