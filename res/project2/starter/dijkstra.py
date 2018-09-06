#### Part 1 functions from PCRS ####

def get_closest(unvisited):
    """ Return a tuple from list unvisited that has the shortest distance
    (return the first such tuple in case of ties).  Assumes unvisited is not
    empty.
    
    Parameters:
        unvisited-- a list of (city_name, distance) pairs
    Return value:
        the first tuple in unvisited whose distance is minimum
    Side-effects:
        None
    
    Examples:
    >>> get_closest([('A', 3)])
    ('A', 3)
    >>> get_closest([('A', 3), ('B', 2)])
    ('B', 2)
    >>> get_closest([('A', 3), ('C', 4)])
    ('A', 3)
    >>> get_closest([('A', 3), ('B', 2), ('C', 4), ('D', 2)])
    ('B', 2)
    """
    
    # Examine every tuple and remember the index of the one with the smallest
    # distance so far-- starting with the first tuple.


def find_city(city, city_list):
    """ (str, list of (str, int)) -> int
    
    Return the index of the first tuple that contains city in city_list, or
    -1 if city does not appear in city_list.
    
    Parameters:
        city-- the name of a city to look for (a string)
        city_list-- a list of tuples of the form (city_name, distance)
    Return value:
        the index of the first tuple in city_list that contains city; -1 if no
        tuple in city_list contains city
    Side-effects:
        None
    
    Examples:
    >>> find_city('A', [('A', 2)])
    0
    >>> find_city('A', [('B', 3), ('C', 2)])
    -1
    >>> find_city('A', [('B', 3), ('C', 2), ('A', 2)])
    2
    >>> find_city('C', [('B', 3), ('C', 2), ('A', 2), ('C', 4)])
    1
    """


def process_line(line):
    """ (str) -> (str, str, int) 
    
    Process one line from a data file (in the same format as the lines in cities.txt)
    to extract the first city's name, the second city's name, and the distance.  Return
    these values in a tuple.
    
    Parameters:
        line-- a single line in the same format as the lines in cities.txt, in the 
               format specified in the project handout (first city:second city distance)
    Return value:
        (first, second, distance), where first is the name of the first city,
        second is the name of the second city, and distance is the distance
    Side-effects:
        None
    
    Examples:
    >>> process_line("A:B 3")
    ('A', 'B', 3)
    >>> process_line("some city:B 2")
    ('some city', 'B', 2)
    >>> process_line("A:other city 5")
    ('A', 'other city', 5)
    """

    
def build_adjacent_distances(lines):
    """ Read distances between cities from the lines read from data,
    which is in the same format as cities.txt, and
    return a dictionary structure that contains all of this information.
    
    Parameters:
        lines -- a list of lines read in from data in the same format as
                 cities.txt
    Return value:
        a dictionary whose keys are city names and whose values are lists of
        pairs of the form (city_name, distance).
        The cities must be sorted in alphabetical order. 
        (Note: sorted([('B', 3), ('A', 2)]) returns [('A', 2), ('B', 3)]
    
    Side-effects:
        None
    
    Examples:
    # The second line below violates style guidelines (it is too long), but
    # is shown this way so you can use it for testing.
    >>> build_adjacent_distances(lines)  # assuming lines are the lines of cities.txt
    {'Toronto': [('New York', 3), ('Mexico City', 7), ('San Francisco', 6)], 'San Francisco': [('Washington', 5), ('Mexico City', 3), ('Toronto', 6)], 'New York': [('Toronto', 3), ('Washington', 2)], 'Washington': [('New York', 2), ('San Francisco', 5)], 'Mexico City': [('San Francisco', 3), ('Toronto', 7)]}
    """


def get_all_cities(city_to_city_dist):
    """ Return a list of all the cities that appear in the dictionary
    city_to_city_dist (a dictionary in the format returned by build_adjacent_distances).
    The cities are to be sorted in alphabetic order.
    
    >>> city_to_city_dist = {'Washington': [('New York', 2), ('San Francisco', 5)],
                    'Toronto': [('New York', 3)]}
    >>> get_all_cities(city_to_city_dist)
    ['New York', 'San Francisco', 'Toronto', 'Washington'] 
    """

    
#### Similar to Part 2 functions on PCRS ####
    
def visit_next(visited, unvisited, city_to_city_dist):
    """ Move the next closest city from the unvisited list to the visited list.
    Update the distances in the unvisited list if a shorter path exists to one
    of the unvisited cities, as described in the handout.  Assumes that the
    unvisited list is non-empty.
    
    Parameters:
        visited-- a list of tuples for cities that have been "visited", i.e.,
            their minimum distance to the city of origin is known, and their
            neighbours already belong to the visited or unvisited list
        unvisited-- a list of tuples for cities that have not yet been visited
        city_to_city_dist-- a dictionary of direct flight lengths between cities, 
                    such as what is returned by get_all_cities()
    Return value:
        None
    Side-effects:
        - the first city C whose distance is minimum in unvisited is removed
          from unvisited and added to visited
        - neighbours of C that did not already belong to either list are added
          to unvisited (with their distance from C)
        - neighbours of C that were already in unvisited have their distance
          updated, if going through C leads to a shorter total distance
    
    Examples:
    # This needs to be tested much more thoroughly than the few cases below, to
    # take into account all the possible situations that could come up.
    >>> city_to_city_dist = build_adjacent_distances(open("cities.txt").readlines())
    >>> unvisited = [('Toronto', 0)]
    >>> visited = []
    >>> visit_next(visited, unvisited, city_to_city_dist)
    >>> visited
    [('Toronto', 0)]
    >>> unvisited
    [('New York', 3), ('Mexico City', 7), ('San Francisco', 6)]
    >>> visit_next(visited, unvisited, city_to_city_dist)
    >>> visited
    [('Toronto', 0), ('New York', 3)]
    >>> unvisited
    [('Mexico City', 7), ('San Francisco', 6), ('Washington', 5)]
    >>> visit_next(visited, unvisited, city_to_city_dist)
    >>> visited
    [('Toronto', 0), ('New York', 3), ('Washington', 5)]
    >>> unvisited
    [('Mexico City', 7), ('San Francisco', 6)]
    """
    
    # Step 1:
    # Find the closest city in the unvisited list and move it to the visited list.
    
    ############################################################################
    
    # Step 2:
    # For every neighbour of the city that was moved to visited, there are three
    # cases to consider:
        # ...(a) neighbour is already visited: do nothing
        # ...(b) neighbour is already unvisited: update its distance
        # ...(c) neighbour is not even unvisited: add it to unvisited
                

def visit_all(city, city_to_city_dist):
    """Return the list of shortest distances from city to every city in 
    the dictionary city_to_city_dist. This is accomplished by repeatedly calling
    visit_next inside a while loop.
    
    
    Parameters:
            city--the starting point
            city_to_city_dist-- a dictionary of direct flight lengths between cities, 
            such as what is returned by get_all_cities()
    
    >> city_to_city_dist = build_adjacent_distances(open("cities.txt").readlines())
    >> visit_all('Toronto', city_to_city_dist)
    [('Mexico City', 7),
     ('New York', 3),
     ('San Francisco', 6),
     ('Toronto', 0),
     ('Washington', 5)]
    """

    
    

    
    
    

            
    
    
