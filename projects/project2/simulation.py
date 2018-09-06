#### Functions from Part 3 on PCRS ####

def get_cities(prob_pairs):
    """Return the list of cities that appear in the list of pairs prob_pairs,
    in the same order that they appear in the list of pairs prob_pairs
    
    Arguments:
      prob_pairs -- a list of city-distance pairs
      
    >> prob_pairs = [('New York', 0),
                     ('Washington', 2),
                     ('Toronto', 3),
                     ('San Francisco', 7),
                     ('Mexico City', 10)]
    >> get_cities(prob_pairs)
    ['New York', 'Washington', 'Toronto', 'San Francisco', 'Mexico City']
    
    """


def get_probabilities(prob_pairs):
    """Return the list of probabilities that appear in the list of pairs
    prob_pairs, in the same order that they appear in prob_pairs.
    
    Arguments:
      prob_pairs -- a list of city-distance pairs
      
    >> prob_pairs = [('New York', 0.1),
                     ('Washington', 0.2),
                     ('Toronto', 0.3),
                     ('San Francisco', 0.2),
                     ('Mexico City', 0.2)]
    >> get_probs(prob_pairs)
    [0.1, 0.2, 0.3, 0.2, 0.2]
    
    """


def init_zero_sick_population(cities):
    """Return a dictionary whose keys are the values in the list cities,
    and whose values are all 0
    
    >> init_zero_sick_population(['TO', 'NYC'])
    {'TO': 0, 'NYC': 0}
    
    Arguments:
      cities -- a list of strings
    
    """



def build_transition_probs(city_to_city_dist, alpha):
    """Return a dictionary representing the probabilities of moving
    from all cities to all other cities according to the formula in
    the project handout.
    
    Arguments: 
      all_dists -- a city_to_city_dist whose keys are cities, and whose
                   values are lists of city-distance pairs
      alpha     -- a float
      
    >> city_to_city_dist = {'Mexico City': [('Mexico City', 0),
                    ('San Francisco', 3),
                    ('Toronto', 7),
                    ('Washington', 8),
                    ('New York', 10)],
                    'New York': [('New York', 0),
                    ('Washington', 2),
                    ('Toronto', 3),
                    ('San Francisco', 7),
                    ('Mexico City', 10)],
                    'Toronto': [('Toronto', 0),
                    ('New York', 3),
                    ('Washington', 5),
                    ('San Francisco', 6),
                    ('Mexico City', 7)],
                    'San Francisco': [('San Francisco', 0),
                    ('Mexico City', 3),
                    ('Washington', 5),
                    ('Toronto', 6),
                    ('New York', 7)],
                    'Washington': [('Washington', 0),
                    ('New York', 2),
                    ('San Francisco', 5),
                    ('Toronto', 5),
                    ('Mexico City', 8)]}

    >> alpha = 2
    
    >> build_transition_probs(city_to_city_dist, alpha)
                {'Mexico City': [('Mexico City', 0.9101374498147844),
                ('San Francisco', 0.056883590613424025),
                ('Toronto', 0.014220897653356006),
                ('Washington', 0.011236264812528202),
                ('New York', 0.007521797105907309)],
                'New York': [('New York', 0.8350726686715951),
                ('Washington', 0.09278585207462167),
                ('Toronto', 0.052192041791974696),
                ('San Francisco', 0.013048010447993674),
                ('Mexico City', 0.0069014270138148355)],
                'Toronto': [('Toronto', 0.8878542892195415),
                ('New York', 0.05549089307622134),
                ('Washington', 0.02466261914498726),
                ('San Francisco', 0.018119475290194722),
                ('Mexico City', 0.013872723269055335)],
                'San Francisco': [('San Francisco', 0.8878542892195415),
                ('Mexico City', 0.05549089307622134),
                ('Washington', 0.02466261914498726),
                ('Toronto', 0.018119475290194722),
                ('New York', 0.013872723269055335)],
                'Washington': [('Washington', 0.8481675392670158),
                ('New York', 0.09424083769633508),
                ('San Francisco', 0.02356020942408377),
                ('Toronto', 0.02356020942408377),
                ('Mexico City', 0.010471204188481676)]}

    """
    

#### Part 4 ####

def build_shortest_distances(direct_distance_dict):    
    """Return a dictionary all_dists whose keys are cities which appear in 
    direct_distance_dict and whose values are the shortest path distance to
    every other city (i.e., the values are return values of visit_all).
    
    
    Arguments:
            direct_flight_distances-- a dictionary of direct flight lengths between cities, 
            such as what is returned by build_adjacent_distances()
    
    
    >> distances = build_adjacent_distances(open("cities.txt").readlines())
    >> build_shortest_distances(distances)
        {'Mexico City': [('Mexico City', 0),
        ('San Francisco', 3),
        ('Toronto', 7),
        ('Washington', 8),
        ('New York', 10)],
        'New York': [('New York', 0),
        ('Washington', 2),
        ('Toronto', 3),
        ('San Francisco', 7),
        ('Mexico City', 10)],
        'Toronto': [('Toronto', 0),
        ('New York', 3),
        ('Washington', 5),
        ('San Francisco', 6),
        ('Mexico City', 7)],
        'San Francisco': [('San Francisco', 0),
        ('Mexico City', 3),
        ('Washington', 5),
        ('Toronto', 6),
        ('New York', 7)],
        'Washington': [('Washington', 0),
        ('New York', 2),
        ('San Francisco', 5),
        ('Toronto', 5),
        ('Mexico City', 8)]}
 
    """

    
def time_step(sick_pop, transition_probs, beta, gamma):
    """
    Change dictionary sick_pop to account for one time-step in the simulation.
    
    Arguments:
    sick_pop: dictionary of city to num sick
    transition_probs: dictionary city to list of (city, prob of destination)
    beta: recovery probability
    gamma: infection probability
    """

###############################################################################
# Sample run of the simulation    
    
alpha = 2       #larger alpha => smaller likelihood of transition
beta = .3       #recovery probability
gamma = .3      #infection probability

# load the dictionary of distances of direct flights between cities
distances = build_adjacent_distances(open("cities.txt").readlines())

# build a dictionary of the shortest distances from every city to every other
all_dists = build_shortest_distances(distances)

# build a dictionary of probabilities of moving from every city to every other
all_probs = build_transition_probs(all_dists, alpha)

# get a list of all the cities in our simulation
cities = get_cities(distances)

# sick_pop is a dictionary of cityname to number of sick inhabitants
# initially every city has 0 sick people
sick_pop = init_zero_sick_population(cities)
# but Toronto has 1 sick person
sick_pop["Toronto"] = 1

# now run the simulation for 100 time steps
for i in range(100):
    time_step(sick_pop, all_probs, beta, gamma)
    print("Day", i, sick_pop)
    

