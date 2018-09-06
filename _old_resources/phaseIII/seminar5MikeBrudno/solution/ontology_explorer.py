import ontology_parser as op

def get_single_path(id_to_parents, start, end):
    """ (dict of {str: list of str}, str, str) -> list of str

    Return a list of the IDs for a single path from start to end
    in the ontology represented by id_to_parents.

    For example, for start '0100271' and end '0000118', return:
    ['0100271', '0001608', '0000118']

    Because a phenotypic abnormality can have multiple parents, there
    can be multiple paths from start to end.  In that case, return any one path.
    """

    # TODO begin
    curr = start
    path = [curr]
    #print(curr, id_to_name[curr])
    while curr != end:
        next_nodes = id_to_parents[curr]
        curr = next_nodes[0]
        #print(curr, id_to_name[curr])
        path.append(curr)
    return path
    # TODO end

## OMIT FROM STARTER begin
def invert_dictionary(id_num_to_parents):
    parent_to_children = {}
    for id_num in id_num_to_parents:
        parents = id_num_to_parents[id_num]
        for parent in parents:
            if parent not in parent_to_children:
                parent_to_children[parent] = [id_num]
            else:
                parent_to_children[parent].append(id_num)
    return parent_to_children 
## OMIT FROM STARTER end

def get_ids_at_level(id_to_parents, level_num, root_id):
    """ (dict of {str: list of str}, str, int) -> list of str

    Return a list of the IDs from level level_num of the
    ontology represented by id_to_parents with root root_id.
    Each ID should appear in the list at most once.

    For example, for level_num  1 and root_id '0000118',
    this function should return:
    ['0001626', '0025031', '0001871', '0001939', '0000924', '0001608', 
    '0040064', '0000707', '0002664', '0001507', '0000769', '0000152', 
    '0003549', '0000478', '0001197', '0000598', '0045027', '0001574', 
    '0000818', '0002715', '0000119', '0002086', '0003011']
    """

    ### TODO begin
    parent_to_children = invert_dictionary(id_to_parents)

    curr = 0
    curr_nodes = [root_id]
    while curr != level_num:
        next_nodes = []
        for id_num in curr_nodes:
            if id_num in parent_to_children:
                next_nodes.extend(parent_to_children[id_num])
        curr += 1 
        curr_nodes = set(next_nodes)
        #print('LEVEL', curr)
        #print(curr_nodes)
    return list(curr_nodes)
    ### TODO end


def get_all_paths(id_to_parents, start, end, path=[]):
    """ (dict of {str: list of str}, str, int) -> list of str

    Return list of list of str representing paths from start to end
    in id_to_parents, including the IDs start and end.

    For example, id_num 0002813 and root_id 0000118, return:
    [['0002813', '0011844', '0011842', '0000924', '0000118'], 
     ['0002813', '0040068', '0000924', '0000118'], 
     ['0002813', '0040068', '0040064', '0000118']]
    """

    # Add start to the current path.
    path = path + [start]

    # Base case: If start and end are equal, return [[start]].
    if start == end:
        return [path]

    # Base case: If start does not have any parents, return the empty list.
    if not start in id_to_parents:
        return []

    # A list of paths (a list of list of str).
    paths = []

    # Recursive case: Use recursion to populate paths.
    for node in id_to_parents[start]:
        if node not in path:
            newpaths = get_all_paths(id_to_parents, node, end, path)
            for newpath in newpaths:
                paths.append(newpath)

    # Return a list of paths (a list of list of str).
    return paths


    # TODO begin
    # if id_num == root_id:
    #     print(id_num, ontology.get_name(id_num))
    #     print()
    # else:
    #     for parent in ontology.get_parents(id_num):
    #         print(id_num, ontology.get_name(id_num))
    #         print_all_paths_to_root(ontology, parent, root_id)
    # TODO end


def get_all_paths_to_level(id_to_parents, start, level_ids, path=[]):
    """ (dict of {str: list of str}, str, int) -> list of str

    Return list of list of str representing paths from start to 
    the level containing level_ids ids in id_to_parents.

    For example, id_num 0009122 and root_id 0000118, and
    the level ids from level 2, return:
    [['0009122', '0009115', '0011842'], ['0009122', '0009121', '0011842']]
    """
    path = path + [start]
    if start in level_ids:
        return [path]
    if not start in id_to_parents:
        return []
    paths = []
    for node in id_to_parents[start]:
        if node not in path:
            newpaths = get_all_paths_to_level(id_to_parents, node, level_ids, path)
            for newpath in newpaths:
                paths.append(newpath)
    return paths
            

if __name__ == '__main__':

    # Read and parse human phenotips ontology.
    ontology_file = open('hp.obo', encoding='utf-8')
    ontology = op.Ontology(ontology_file)

    # Get the ontology dictionaries.
    id_to_parents = ontology.get_id_to_parents()
    id_to_name = ontology.get_id_to_name()

    # ID of the root of the phenotypic abnormability subcategory.
    root = '0000118'
    
    # == Example calls on the required functions == 

    # For id_num '0100271', there is only one path to root.
    #print(get_single_path(id_to_parents, '0100271', root))

    # For id_num '0009122', there is are multiple paths to root,
    # since it has two parents.  Print one of them.
    #print(get_single_path(id_to_parents, '0009122', root))
    
    # Phenotypic abnormalities by level (distance from root).
    #print(get_ids_at_level(id_to_parents, 1, root))
    #print(get_ids_at_level(id_to_parents, 3, root))

    # Get all paths to root.
    #print(get_all_paths(id_to_parents, '0009122', root))
    #print(get_all_paths(id_to_parents, '0004855', root))

    # Print all paths from 0009122 to level 2.
    #level_ids = get_ids_at_level(id_to_parents, 2, root)
    #print(get_all_paths_to_level(id_to_parents, '0002813', level_ids))
    #print(get_all_paths_to_level(id_to_parents, '0009122', level_ids))

