import functools
import operator

def bitvector(iterable):
    return functools.reduce(operator.or_, (1 << x for x in iterable))
