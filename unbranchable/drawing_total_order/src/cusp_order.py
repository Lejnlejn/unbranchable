# cusp_order.py
# Code for working with edge orientations on custom continents. In particular, we compute total orders on cusps
# restricted to the vertices of a finite veering continent. 

from collections import deque
from .edge_orientations import orient_edges

# We give a mathematical guarantee that would be a distraction from the paper, so we state it here.
#
# Let V be a veering triangulation and L the corresponding loom space.
# We borrow terminology from Schleimer-Segerman 2021 and Parlak-Segerman 2026.
### Lemma: Let C\subset lift(V) be a convex continent. Let Loom(C)\subset L be its loomification. Then if R\subset Loom(C) is a skeletal rectangle, the 
### corresponding cell of lift(V) is contained in C. 
### Proof: Following Definition 7.5 Parlak-Segerman claim (the argument is sufficiently brief): C is convex if and only if C contains all tetrahedron
### whose corresponding rectangles are contained in Loom(C). Any edge rectangle contained in Loom(C) belongs
### to some tetrahedron rectangle by Proposition 6.46 of Schleimer-Segerman. Namely, if this were not the case, then the content of Loom(C) would not
### be a triangulated ball. Since C contains all tetrahedra whose rectangles are in Loom(C), it then contains all edges whose edge rectangles
### are in Loom(C). //
#
# Thank you to Anna Parlak for a suggestion that made the above proof much quicker, and both Anna and Henry for their proofs of the equivalence
# of the definitions.
#
# Let e_i be a path given by edges of lift(V).
# Call e_i _geodesic_ if the path obtained by following a (choose any) material boundary side of each of the edges is a geodesic in L.
#
### Claim: Let C be a convex veering continent. Then for every pair of vertices u, v of C, there is an edge path from u to v, consisting
### of edges in C, that is geodesic.
### Proof sketch: There exists some path between the endpoints by edges. If the path is not geodesic, then there is a U-turn in the path. 
#--------- A U-turn in such a path is a pair of sequential edges e_i and e_{i+1} such that e_i and e_{i+1} enter and exit from the
#--------- same or adjacent orthant of a vertex v. Say that a is the tail of e_i and c the head of e_{i+1}. Then the geodesic hull between a
#--------- and c is a single convex polygon P. There is then a collection of edges d_i of lift(V) from u to w whose loomification is contained in P.
#--------- These edge rectangles are thuse contained in Loom(C). By Lemma, C contains these edges. 
#--------- 
#--------- Now we replace the subset of the sequence e_i and e_{i+1} with the sequence d_i. There are no U-turns contained in this d_i subset, 
#--------- so we decrease a complexity in a manner analogous to Lemma 6.11 of the SS's loom space paper. //
#
### With above claim, we can always find geodesics in a convex continent. Note that one does need the above claim to "try" to find the geodesics;
### the after-the-fact computation is valid regardless. 

def compute_geodesic_data(ball):
    # Builds all the data used to compute whether or not an edge path is geodesic.

    if ball.geodesic_data is not None: # if we have already computed, move on
        return ball.geodesic_data

    vertex_pair_to_edge = {} # (v0,v1)->edge num of edge between them, if it exists 
    adjacency = {vertex: set() for vertex in range(ball.tri.countVertices())} # the underlying adjacency relations given by the edges
    incident = {vertex: [] for vertex in range(ball.tri.countVertices())}

    for edge in ball.tri.edges():
        vertex_a = edge.vertex(0).index()
        vertex_b = edge.vertex(1).index()
        if vertex_a == vertex_b:
            continue
        edge_index = edge.index()
        vertex_pair_to_edge[frozenset((vertex_a, vertex_b))] = edge_index
        adjacency[vertex_a].add(vertex_b)
        adjacency[vertex_b].add(vertex_a)
        incident[vertex_a].append((vertex_b, edge_index))
        incident[vertex_b].append((vertex_a, edge_index))

    positions = {vertex: index for index, vertex in enumerate(ball.vert_circle)}
    circle_size = len(ball.vert_circle)

    
    colour_fan_of_edge_at_vertex = {} # (vertex,edge)-> which numbered colour fan at the edge of verte xedge lies in

    for vertex in range(ball.tri.countVertices()):
        if vertex not in positions:
            continue
        
        vertex_position = positions[vertex]
        # the fan is ordered by the position of the other vertex of the edge on the circle cut at vertex
        fan = sorted(incident[vertex], key=lambda pair: (positions[pair[0]] - vertex_position) % circle_size)
        colour_runs = []
        for other_vertex, edge_index in fan:
            # continue the same colour fan
            colour = ball.edge_colour(edge_index)
            if colour_runs and colour_runs[-1][1] == colour: 
                colour_runs[-1][0].append(edge_index)
            else:
                # make a new fan for the altered colour
                colour_runs.append(([edge_index], colour))

        for run_index, (run_edges, run_colour) in enumerate(colour_runs):
            for edge_index in run_edges:
                colour_fan_of_edge_at_vertex[vertex, edge_index] = run_index

    data = {'vertex_pair_to_edge': vertex_pair_to_edge, 'adjacency': adjacency, 'colour_fan_of_edge_at_vertex': colour_fan_of_edge_at_vertex}
    ball.geodesic_data = data
    return data

def turn_is_geodesic(vertex, edge_before, edge_after, colour_fan_of_edge_at_vertex):
    run_before = colour_fan_of_edge_at_vertex.get((vertex, edge_before))
    run_after = colour_fan_of_edge_at_vertex.get((vertex, edge_after))
    return run_before is not None and run_after is not None and abs(run_before - run_after) >= 2

# Observe that a geodesic edge path is equivalent to the following: at every cusp, the orthants we enter and leave from
# are distinct and non-adjacent. Moreover, note that the orthants are deliniated by continuous fans of the same colour.

def find_geodesic_path(ball, source, target):
    # Find a geodesic path using BFS
    data = compute_geodesic_data(ball)
    adjacency = data['adjacency']
    colour_fan_of_edge_at_vertex = data['colour_fan_of_edge_at_vertex']
    vertex_pair_to_edge = data['vertex_pair_to_edge']

    starting_path = [source]
    starting_last_edge = None
    starting_visited_vertices = frozenset({source})

    queue = deque([(starting_path, starting_last_edge, starting_visited_vertices)])

    while queue:
        vertex_path, last_edge, visited_vertices = queue.popleft()
        current_vertex = vertex_path[-1]

        if current_vertex == target:
            return vertex_path
        
        # explore to neighbours
        for next_vertex in adjacency[current_vertex]:
            
            if next_vertex in visited_vertices:
                continue

            next_edge = vertex_pair_to_edge[frozenset((current_vertex, next_vertex))]

            if last_edge is not None:
                # check if we are continuing along a geodesic
                if not turn_is_geodesic(current_vertex, last_edge, next_edge, colour_fan_of_edge_at_vertex):
                    continue

            new_path = vertex_path + [next_vertex]
            new_visited_vertices = visited_vertices | {next_vertex}

            queue.append((new_path, next_edge, new_visited_vertices))
            
    assert False # no geodesic found

#### Now we turn to computing geodesic windings with the geodesic representatives.

def orient_ball_edges(ball):
    # returns a list of edge orientations for the continent ball. 
    # TODO: we should simply lift these from downstairs; this is a continent-sized computation.
    coors = ball.coors2

    if coors is None or any(coor is None for coor in coors):
        coors = None

    return orient_edges(ball.tri, ball.angle2, ball.merge2, tet_vert_coorientations=coors)

def direction_lookup(ball, edge_orientations):
    # pairs of vertices (v0,v1) as edge endpoints -> if the edge orientation points with that ordered tuple forward
    forwards = {}

    for edge in ball.tri.edges():
        edge_index = edge.index()
        vertex_a = edge.vertex(0).index()
        vertex_b = edge.vertex(1).index()

        if vertex_a != vertex_b:

            if edge_orientations[edge_index] == 1:
                forwards[vertex_a, vertex_b] = True
                forwards[vertex_b, vertex_a] = False
            else:
                forwards[vertex_a, vertex_b] = False
                forwards[vertex_b, vertex_a] = True

    return forwards

def subgeodesics_from_vertex_path(ball, edge_orientations, vertex_path):
    # converts a vertex_path into "subgeodesics": a list of (vertex_from, vertex_to, edge_num, are we following the edge forwards)
    data = compute_geodesic_data(ball)
    forwards = direction_lookup(ball, edge_orientations)
    subgeodesics = []
    for vertex_from, vertex_to in zip(vertex_path[:-1], vertex_path[1:]):
        edge_index = data['vertex_pair_to_edge'][frozenset((vertex_from, vertex_to))]
        subgeodesics.append((vertex_from, vertex_to, edge_index, forwards[vertex_from, vertex_to]))

    return subgeodesics

# Now we want to find paths of edges between two vertices.
def chain_directions(subgeodesics):
    # counts the signs in the subgeodesic decomposition of the edge-path geodesic
    directions = []
    previous = None
    for subgeodesic in subgeodesics:
        forwards = subgeodesic[3]
        if forwards != previous:
            directions.append(forwards)
            previous = forwards

    return directions

jump_assigner = {'left': -1, 'right': 1}

def twisted_geodesic_windings(ball, subgeodesics):
    # Jump signs along a path.
    # This is called the untwisted geodesic winding at the vertex in the paper.
    signs = []
    for index in range(len(subgeodesics) - 1):
        vertex_in, vertex_jump, edge_in, forwards_in = subgeodesics[index]
        vertex_jump_out, vertex_out, edge_out, forwards_out = subgeodesics[index + 1]

        if forwards_in != forwards_out:
            assert vertex_jump == vertex_jump_out  # consecutive subgeodesics meet
            planar = jump_assigner[ball.left_or_right((vertex_in, vertex_jump), vertex_out)]
            incoming = 1 if forwards_in else -1
            signs.append(planar * incoming) # here is the twisted modification, * incoming, from untwisted. 

    return signs

def score_subgeodesics(ball, subgeodesics):
    # some of signs along subgeodesic divsion. This is the "sigma" term. 
    chain_sum = sum(-1 if forwards else 1 for forwards in chain_directions(subgeodesics))
    return chain_sum + sum(twisted_geodesic_windings(ball, subgeodesics))

# Sort all vertices from scratch using geodesic comparisons.
def Theta(ball, edge_orientations, comparison_data, vertex_a, vertex_b):
    # Apply the order function from the paper to the pair vertex_a and vertex_b
    if (vertex_a, vertex_b) in comparison_data:
        return comparison_data[vertex_a, vertex_b]
    
    vertex_path = find_geodesic_path(ball, vertex_a, vertex_b)
    assert vertex_path # no geodesic found
    score = score_subgeodesics(ball, subgeodesics_from_vertex_path(ball, edge_orientations, vertex_path))
    comparison_data[vertex_a, vertex_b]= score
    comparison_data[vertex_b, vertex_a]= -score

    return score

def compare(ball, edge_orientations, comparison_data, vertex_a, vertex_b):
    # Compute the verdict of the total order according to Theta
    score = Theta(ball, edge_orientations, comparison_data, vertex_a, vertex_b)
    return 1 if score > 0 else -1 if score < 0 else 0

def winding_sign(ball, edge_orientations, vertex_a, vertex_b):
    # Extend a previous order by binary insertion.
    vertex_path = find_geodesic_path(ball, vertex_a, vertex_b)
    assert vertex_path is not None # no geodesic found
    score = score_subgeodesics(ball, subgeodesics_from_vertex_path(ball, edge_orientations, vertex_path))

    return 1 if score > 0 else -1 if score < 0 else 0

def update_order_incremental(ball, edge_orientations, previous_order=None, old_vertex_count=None):
    # this is the slow way of computing the total order. It literally runs the geodesic comparision on every pair.
    compute_geodesic_data(ball)
    vertex_count = ball.tri.countVertices()

    if previous_order is None:
        scores, ranks = compute_scores_via_sort(ball, edge_orientations)
        return scores, ranks, sorted(range(vertex_count), key=lambda vertex: ranks[vertex])

    order = list(previous_order)
    comparisons = 0
    for vertex in range(old_vertex_count, vertex_count):
        # TODO: add an exit clause if this fails; it means that we do not have a total order. 
        low = 0
        high = len(order)
        while low < high:
            middle = (low + high) // 2
            comparisons += 1
            if winding_sign(ball, edge_orientations, vertex, order[middle]) > 0:
                low = middle + 1
            else:
                high = middle
        order.insert(low, vertex)
    ranks = {vertex: index for index, vertex in enumerate(order)}
    scores = {(vertex_a, vertex_b): 1 if ranks[vertex_a] > ranks[vertex_b] else -1 for vertex_a in range(vertex_count) for vertex_b in range(vertex_count) if vertex_a != vertex_b}

    return scores, ranks, order

# Given our proof we know that Theta defines a total order, we can shortcut the total computation by taking closures of
# the relations we already know.
def compute_scores_via_sort(ball, edge_orientations):
    compute_geodesic_data(ball)
    vertex_count = ball.tri.countVertices()
    comparison_data = {}

    # Seed the comparisons using edges; these are quick consequences
    for edge in ball.tri.edges():
        vertex_a = edge.vertex(0).index()
        vertex_b = edge.vertex(1).index()

        if vertex_a == vertex_b or (vertex_a, vertex_b) in comparison_data:
            continue

        subgeodesics = subgeodesics_from_vertex_path(ball, edge_orientations, [vertex_a, vertex_b])
        score = score_subgeodesics(ball, subgeodesics)

        comparison_data[vertex_a, vertex_b] = score
        comparison_data[vertex_b, vertex_a] = -score

    # Insert each vertex into its correct place in the total order
    order = []
    for vertex in range(vertex_count):
        # following on the current order until we find our spot. The closure relation allows us to just check locally.
        position = 0

        while position < len(order) and compare(ball, edge_orientations, comparison_data, vertex, order[position]) > 0:
            position += 1

        order.insert(position, vertex)
        # print(order)

    ranks = {vertex: index for index, vertex in enumerate(order)}
    scores = {(vertex_a, vertex_b):1 if ranks[vertex_a] > ranks[vertex_b] else -1 for vertex_a in range(vertex_count) for vertex_b in range(vertex_count) if vertex_a != vertex_b}

    return scores, ranks