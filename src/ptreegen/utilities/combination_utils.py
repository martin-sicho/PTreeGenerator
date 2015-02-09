def uniqueCombinationsGenerator(items, n):
    if not n:
        yield []
    else:
        for i in xrange(len(items)):
            for comb in uniqueCombinationsGenerator(items[i+1:], n-1):
                yield [items[i]] + comb

def combinationsGenerator(items, n):
    if not n:
        yield []
    else:
        for i in xrange(len(items)):
            for comb in combinationsGenerator(items[:i]+items[i+1:],n-1):
                yield [items[i]] + comb