## @package combination_utils
# Contains functions that generate
# different kinds of combinations
# from an iterable.
#

##
# A generator of unique combinations of size n
# from the given list of items.
#
# @param items an iterable to be processed
# @param n size of each of the unique combinations generated
# @return returns all possible unique combinations as a list of lists
def uniqueCombinationsGenerator(items, n):
    if not n:
        yield []
    else:
        for i in xrange(len(items)):
            for comb in uniqueCombinationsGenerator(items[i+1:], n-1):
                yield [items[i]] + comb

##
# A generator of all non-unique combinations
# (permutations) of size n from the given list of items.
#
# @param items an iterable to be processed
# @param n size of each of the combinations generated
# @return returns all possible unique combinations as a list of lists
def combinationsGenerator(items, n):
    if not n:
        yield []
    else:
        for i in xrange(len(items)):
            for comb in combinationsGenerator(items[:i]+items[i+1:],n-1):
                yield [items[i]] + comb