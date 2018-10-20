def search_v1(lst, key):
    for item in lst:
        if item == key:
            return True
    return False

def search_v2(lst, key): 
    if lst == []:       # base case: list is empty 
        return False
    elif lst[0] == key: # base case: found key
        return True
    else:               # recursive case: search remainder of list
        return search_v2(lst[1:], key)
