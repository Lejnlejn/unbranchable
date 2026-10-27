# Track data for the landscapes of a custom continent
# This currently assumes that the continent we input is convex, and that we go looking for tracks on the appropriate
# upper/lower boundary of the continent.

from .branchable import face_num_to_edges
from veering.taut import edge_num_to_vert_pair, vert_pair_to_edge_num


# throughout, parity = 'up' or 'down' for upper and lower track, respectively
def large_track_edge_in_face(ball, face_index, parity):
    # gets the large edge of the upper/lower track in the face

    face = ball.tri.triangle(face_index)
    embeddings = list(face.embeddings())
    assert len(embeddings)==1 # we currently only use this on boundary faces
    
    embedding = embeddings[0]
    tetrahedron = embedding.simplex()
    face_num = embedding.face()
    if parity=='up':
        target_coor = 1
    else:
        target_coor = -1

    if ball.coors[tetrahedron.index()][face_num] != target_coor:
        assert False # looking for the wrong track on the wrong boundary
    
    downstairs_tetrahedron = ball.lifteds[tetrahedron.index()].home
    downstairs_face = ball.base_tri.triangle(ball.face_to_quotient(face_index))
    # find the other face downstairs
    matching = next(face_embedding for face_embedding in downstairs_face.embeddings() if face_embedding.simplex() == downstairs_tetrahedron and face_embedding.face() == face_num)
    other = next(face_embedding for face_embedding in downstairs_face.embeddings() if face_embedding.simplex() != downstairs_tetrahedron or face_embedding.face() != face_num)
    track_tetrahedron = other.simplex()
    track_face_num = other.face()
    pi_edges = (ball.base_angle[track_tetrahedron.index()], 5 - ball.base_angle[track_tetrahedron.index()])
    large_edge = next(pi_edge for pi_edge in pi_edges if pi_edge in face_num_to_edges(track_face_num))
    large_vertices = edge_num_to_vert_pair[large_edge]
    # now map it back over to our face
    gluing = matching.vertices() * other.vertices().inverse()
    transferred_vertices = tuple(gluing[vertex] for vertex in large_vertices)
    transferred_edge = vert_pair_to_edge_num[tuple(sorted(transferred_vertices))]
    return tetrahedron.edge(transferred_edge).index()

class LandscapeTrack:
    def __init__(self, ball, parity):
        assert parity in ('up', 'down')  # upper or lower track
        self.parity = parity
        self.coastal = set(ball.find_coastal_edges()) # the vertices at coastal edges
        self.vertices = set()
        self.edges = []
        for face_index in ball.get_upper_faces(parity=parity):
            large_edge = large_track_edge_in_face(ball, face_index, parity)
            face = ball.tri.triangle(face_index)
            face_edges = [face.edge(index).index() for index in range(3)]
            small_edges = [edge_index for edge_index in face_edges if edge_index != large_edge]
            for small_edge in small_edges:
                self.vertices.add(small_edge)
                self.vertices.add(large_edge)
                self.edges.append((small_edge, large_edge))

    def coastal_track_arcs(self):
        # Coastal endpoints of the long track arcs
        outgoing = {} # small_edge->large_edge options. For internal sinks this will be 2 edges, for everybody else, it is one.  
        valence = {vertex: 0 for vertex in self.vertices}
        for small_edge, large_edge in self.edges:
            if small_edge not in outgoing:
                outgoing[small_edge] = []

            outgoing[small_edge].append(large_edge)

            valence[small_edge] += 1
            valence[large_edge] += 1

        def follow_to_coast(first):
            # follows the track to the coast
            current = first
            steps = 0
            while current not in self.coastal:
                if current not in outgoing:
                    return None
                current = outgoing[current][0]
                steps += 1
                assert steps <= len(self.vertices)  # track does not cycle

            return current

        arcs = []
        for source in self.vertices:
            if source in self.coastal and source in outgoing: # coastal sources
                sink = follow_to_coast(outgoing[source][0])
                if sink is not None and sink != source:
                    arcs.append((source, sink))

        for source in self.vertices:
            if source in self.coastal or valence[source] != 2:
                continue
            
            # Otherwise, we have an interior source, and we have two ways to follow it to the coast.
            if source in outgoing:
                targets = outgoing[source]
            else:
                targets = []

            # if len(targets) != 2:
            #     continue
            assert len(targets)==2 # we are in the valence two case, so we should always have two outgoings.
            
            coast_a = follow_to_coast(targets[0])
            coast_b = follow_to_coast(targets[1])

            if coast_a is not None and coast_b is not None and coast_a != coast_b:
                arcs.append((coast_a, coast_b))

        return arcs