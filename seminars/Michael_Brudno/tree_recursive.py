def get_leaves(tree, node, leaves=[]):
    if node not in tree:  # base case
        leaves.append(node)
    else:                 # recursive case
        for child in tree[node]:
            get_leaves(tree, child, leaves)
        return leaves

tree = {'A': ['B', 'C', 'D'],
     'B': ['E', 'F'],
     'C': ['G', 'H'],
     'H': ['J'],
     'D': ['I']}
print(get_leaves(tree, 'A'))

