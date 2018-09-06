import dijkstra 
    
distances = {'Mexico City': [('San Francisco', 3), ('Toronto', 7)],
                        'New York': [('Toronto', 3), ('Washington', 2)],
                        'San Francisco': [('Mexico City', 3), ('Toronto', 6), ('Washington', 5)],
                        'Toronto': [('Mexico City', 7), ('New York', 3), ('San Francisco', 6)],
                        'Washington': [('New York', 2), ('San Francisco', 5)]}    

unvisited = [('Toronto', 0)]
visited = []
dijkstra.visit_next(visited, unvisited, distances)

if visited == [('Toronto', 0)]:
    print("Test 1 passed")
    
if unvisited == [('Mexico City', 7), ('New York', 3), ('San Francisco', 6)]:
    print("Test 2 passed")
    
dijkstra.visit_next(visited, unvisited, distances)
if visited == [('Toronto', 0), ('New York', 3)]:
    print("Test 3 passed")

if unvisited == [('Mexico City', 7), ('San Francisco', 6), ('Washington', 5)]:
    print("Test 4 passed")

dijkstra.visit_next(visited, unvisited, distances)

if visited == [('Toronto', 0), ('New York', 3), ('Washington', 5)]:
    print("Test 5 passed")

if unvisited == [('Mexico City', 7), ('San Francisco', 6)]:
    print("Test 6 passed")

dijkstra.visit_next(visited, unvisited, distances)

if visited == [('Toronto', 0), ('New York', 3), ('Washington', 5), ('San Francisco', 6)]:
    print("Test 7 passed")
    
dijkstra.visit_next(visited, unvisited, distances)
    
if visited ==  [('Toronto', 0), ('New York', 3), ('Washington', 5), ('San Francisco', 6), ('Mexico City', 7)]:
    print("Test 8 passed")
    
###########################################################################################################################
    
if dijkstra.visit_all("Toronto", distances)  ==  [('Mexico City', 7), ('New York', 3), ('San Francisco', 6), ('Toronto', 0), ('Washington', 5)]:
    print("Test 9 passed")
    
if dijkstra.visit_all("New York", distances)  ==  [('Mexico City', 10), ('New York', 0), ('San Francisco', 7), ('Toronto', 3), ('Washington', 2)]:
    print("Test 10 passed")

if dijkstra.visit_all("San Francisco", distances) == [('Mexico City', 3), ('New York', 7), ('San Francisco', 0), ('Toronto', 6), ('Washington', 5)]:
    print("Test 11 passed")
    
if dijkstra.visit_all("Mexico City", distances) == [('Mexico City', 0), ('New York', 10), ('San Francisco', 3), ('Toronto', 7), ('Washington', 8)]:
    print("Test 12 passed")

if dijkstra.visit_all("Washington", distances) == [('Mexico City', 8), ('New York', 2), ('San Francisco', 5), ('Toronto', 5),  ('Washington', 0)]:
    print("Test 13 passed")
    