#!/usr/bin/env python3

def inthash(key):
  return (key * 14695981039346656037) % (1 << 64)


def short_str(l):
    l = list(l)
    if len(l) == 0:
        return '.'
    return ''.join(map(str,l))


if __name__ == "__main__":
	for i in range(100):
		print(i, inthash(inthash(inthash(i))))
