# Iterative approach.
def reverse_v1(L):
    rev = []
    for item in L:
        rev = [item] + rev
    return rev

# Recursive approach.
def reverse_v2(L):
    if len(L) == 0: # base case
        return []
    else:           # recursive case
        return [L[-1]] + reverse_v2(L[:-1])