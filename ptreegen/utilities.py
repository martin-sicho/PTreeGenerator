from sys import stderr

from ete2 import Tree


def findConsensusTree(trees,weights=[],lim=0):
    """
    Returns weighted consensus tree. Uses 50% majority rule.

    It is a modified version of
    Francois-Jose Serra's code (accesible from here:
    https://github.com/fransua/utils/blob/master/pmodeltest/consensus.py).
    """

    if weights == []:
        weights = [1] * len (trees)
    dic = {}
    outgroup_name = trees[0].get_leaf_names()[1]
    tlen = 0
    for (tree, weight) in zip (trees, weights):
        if tlen == 0: tlen = len(tree)
        elif len (tree) != tlen: exit('ERROR: trees with different length')
        outgroup = tree.search_nodes(name=outgroup_name)[0]
        tree.set_outgroup(outgroup)
        dad = outgroup.get_sisters()[0]
        for node in dad.traverse():
            if node.is_root(): continue
            cluster  = ','.join (sorted (node.get_leaf_names()))
            if dic.has_key(cluster):
                dic[cluster] += weight
            else:
                dic[cluster] =  weight

    sorted_nodes = map(lambda x: [x[2], x[1]], sorted (\
        map (lambda x: (len (x.split(',')), x, dic[x]), \
             dic.keys()), reverse = True))
    if lim < sorted (sorted_nodes, reverse=True)[:tlen*2 - 3][-1][0]:
        lim = sorted (sorted_nodes, reverse=True)[:tlen*2 - 3][-1][0]
    sorted_nodes = filter (lambda x: x[0] >= lim, sorted_nodes)
    sorted_nodes = map (lambda x: x[1], sorted_nodes)
    if len (sorted_nodes) > tlen*2 - 3:
        print >> stderr, \
              'WARNING: two nodes with same support, will remove: ' + \
              sorted_nodes[-1]
        sorted_nodes = sorted_nodes[:-1]
    cons_tree = Tree()
    cons_tree.add_child(name=outgroup_name)
    node = cons_tree.add_child(name='NoName')
    node.add_feature('childrens', \
                     set (sorted_nodes.pop(0).split(','))
                     - set([outgroup_name]))
    while len (sorted_nodes) > 0:
        for name in sorted_nodes:
            if not name in sorted_nodes: continue
            for node in cons_tree.traverse(strategy='postorder'):
                if node.is_root(): continue
                if node.name is not 'NoName': continue
                if len (node.childrens & set(name.split(','))) == 0:
                    continue
                # check if ther is better solution in one of the child
                for rest in sorted_nodes:
                    if len (set(rest.split(','))) < \
                       len (set(name.split(','))):
                        continue
                    if len (set(rest.split(',')) & set(name.split(','))) > 0:
                        name = rest
                weight = dic[name]
                children = set(name.split(','))
                if len (children) == 1:
                    node.add_child(name=name)
                else:
                    n = node.add_child(name='NoName')
                    n.add_feature('childrens', children)
                    n.support = weight
                break
            sorted_nodes.pop(sorted_nodes.index(name))
            sister = node.childrens - children
            name = ','.join (sorted ( list (sister)))
            if not name in sorted_nodes:
                continue
            weight = dic[name]
            if len (sister) == 1:
                node.add_child(name=name)
            else:
                n = node.add_child(name='NoName')
                n.add_feature('childrens', sister)
                n.support = weight
            sorted_nodes.pop(sorted_nodes.index(name))
            break
    # return cons_tree.write(format=9)
    return cons_tree