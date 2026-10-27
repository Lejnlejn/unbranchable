# edge_orientations.py
# Edge orientations from coorientable branchings. The branchings are in the form of a z2-merging structure
# as in branchable.py

from .branchable import edge_side_circle as closed_edge_side_circle, face_num_to_large_edge, is_transverse_taut, split_pairs
from veering.taut import vert_pair_to_edge_num

def edge_side_circle(tri):
    # Face-edge pairs around each edge. This is adapted from branchable for the boundary case. Here, the boundary fans
    # are just ordered linearly as opposed to circularly.
    circles = closed_edge_side_circle(tri) # this is an aliased import so that we just need to handle the boundary cases
    for edge_num in range(tri.countEdges()):
        edge = tri.edge(edge_num)
        
        if not edge.isBoundary():
            continue

        embeddings = order_boundary_edge_embeddings(list(edge.embeddings()))
        first = embeddings[0]
        tetrahedron = first.tetrahedron()
        trailing, leading = first.vertices()[2], first.vertices()[3]
        face = tetrahedron.triangle(trailing)
        edge_in_face = tetrahedron.faceMapping(2, trailing).inverse()[leading]
        pairs = [(face.index(), edge_in_face)]

        for embedding in embeddings:
            tetrahedron = embedding.tetrahedron()
            trailing, leading = embedding.vertices()[2], embedding.vertices()[3]
            face = tetrahedron.triangle(leading)
            edge_in_face = tetrahedron.faceMapping(2, leading).inverse()[trailing]
            pairs.append((face.index(), edge_in_face))

        circles[edge_num] = pairs

    return circles

def trailing_leading_faces(embedding):
    tetrahedron = embedding.tetrahedron()
    vertices = embedding.vertices()
    return tetrahedron.triangle(vertices[2]), tetrahedron.triangle(vertices[3])

def order_boundary_edge_embeddings(embeddings):
    # Order a boundary-edge fan from one boundary face to the other.
    start = next(embedding for embedding in embeddings if trailing_leading_faces(embedding)[0].isBoundary())
    ordered = [start]
    used = {id(start)}
    current = start
    while True:
        trailing_face, leading_face = trailing_leading_faces(current)
        # print(trailing_face, leading_face)
        if leading_face.isBoundary():
            break

        next_embedding = next((embedding for embedding in embeddings if id(embedding) not in used and trailing_leading_faces(embedding)[0] == leading_face), None)
        assert next_embedding is not None  # boundary fan is connected
        ordered.append(next_embedding)
        used.add(id(next_embedding))
        current = next_embedding

    return ordered

def edge_sides_splitted(tri, angle):
    # Split each edge fan at pi edges. This already appears in branchable, but we want specifically the with-boundary case here. 
    edge_sides = edge_side_circle(tri)
    splits = split_pairs(tri, angle)
    edge_chunks = []
    for edge_num in range(tri.countEdges()):
        candidates = list(edge_sides[edge_num])
        splitters = splits[edge_num]
        if len(candidates) == 1:
            edge_chunks.append([[candidates[0]]])
            continue

        chunks = []
        current_chunk = []
        for index in range(len(candidates)):
            next_index = (index + 1) % len(candidates)
            pair_a = candidates[index]
            pair_b = candidates[next_index]
            current_chunk.append(pair_a)
            if [pair_a, pair_b] in splitters or [pair_b, pair_a] in splitters:
                chunks.append(current_chunk)
                current_chunk = []

        if not splitters:
            chunks = [edge_sides[edge_num]]
        elif current_chunk:
            chunks = [current_chunk + chunks[0]] + chunks[1:]

        edge_chunks.append(chunks)

    return edge_chunks

def merge_to_coors(tri, angle, merge, edge_chunks):
    # Merge to coors around each edge
    face_large_edges = [face_num_to_large_edge(face_num, tri, angle, merge) for face_num in range(tri.countTriangles())]
    output = []
    for chunks in edge_chunks:
        edge_coors = []
        for chunk in chunks:
            face_num, edge_in_face = chunk[0]
            edge_coors.append(int(face_large_edges[face_num] == edge_in_face))
        output.append(edge_coors)

    return output

def face_coorientation(face, tet_vert_coorientations):
    face_embedding = face.embedding(0)
    tet = face_embedding.simplex()
    face_index_in_tet = face_embedding.vertices()[3]

    coorientation = tet_vert_coorientations[tet.index()][face_index_in_tet]
    does_face_agree_with_tet = face_embedding.vertices().sign()

    return coorientation * does_face_agree_with_tet


# Now the actual edge orientations.
# A face coorientation is +1 or -1 in the veering convention if the circular order of the orientaiton agrees or disagrees
# with the 0-1-2 ordering on the vertices.

# So, in a +1 face, the tail vertex of an edge is given by:
edge_tail = {0: 1, 1: 2, 2: 0}
edge_head = {0: 2, 1: 0, 2: 1}

def orient_edge_from_face(edge_in_face, edge_coor, face_coor):
    # Local direction induced on one face edge.
    # The right hand-rule computation is the following:
    # If the coorientation points into the edge, orient the edge to agree with the edge orientation given by the face coor.
    # Otherwise, disagree. 

    local_tail = edge_tail[edge_in_face]
    local_head = edge_head[edge_in_face]
    
    if not edge_coor:
        # edge coor points out, flip the edge
        local_tail, local_head = local_head, local_tail
        
    if face_coor == -1:
        # default labelling was not the face labelling, flip the edge
        local_tail, local_head = local_head, local_tail

    return local_tail, local_head

def orient_edges(tri, angle, merge, tet_vert_coorientations=None):
    # Orient all edges from the merge structure. We record a +1 if the oriented edge agree's with regina's embedding of the edge.
    # That is, +1 if the edge is oriented from edge.vertex(0) to edge.vertex(1). Otherwise -1
    # Returns a list of +1/-1s indexed by edge num. 
    

    assert tri.isOriented()  # edge orientation needs oriented tetrahedra
    if tet_vert_coorientations is None:
        tet_vert_coorientations = is_transverse_taut(tri, angle, return_type='tet_vert_coorientations')

    edge_chunks = edge_sides_splitted(tri, angle)
    edge_coors = [coors[0] for coors in merge_to_coors(tri, angle, merge, edge_chunks)]
    orientations = []

    for edge_index in range(tri.countEdges()):
        face_index, edge_in_face = edge_chunks[edge_index][0][0]
        face = tri.triangle(face_index)
        face_coor = face_coorientation(face, tet_vert_coorientations)
        edge_coor = edge_coors[edge_index]
        local_tail, local_head = orient_edge_from_face(edge_in_face, edge_coor, face_coor)
        embedding = face.embedding(0)
        vertices = embedding.vertices()
        tetrahedron = embedding.simplex()

        tet_tail = vertices[local_tail]
        tet_head = vertices[local_head]

        local_edge_num = vert_pair_to_edge_num[min(tet_tail, tet_head), max(tet_tail, tet_head)]
        regina_tail = tetrahedron.edgeMapping(local_edge_num)[0]

        orientations.append(1 if tet_tail == regina_tail else -1)

    return orientations