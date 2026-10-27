# Base-triangulation helpers needed to grow and draw the layered ball.

from .branchable import edge_side_circle
from veering.taut import vert_pair_to_edge_num
from veering.transverse_taut import get_tet_top_vert_nums, edge_side_face_collections

def assert_combed(tri, angle):
    sides = edge_side_face_collections(tri, angle)
    for edge_num, (side0, side1) in enumerate(sides):
        # combed is at least two faces on each side:
        assert len(side0) >= 2 and len(side1) >= 2

def embedge_to_tet_num_pair(embedge):
    return (embedge.simplex().index(), embedge.face())

def find_paired_face(tri, pair):
    tet_num, face_num = pair
    tet = tri.simplex(tet_num)
    adj_tet = tet.adjacentSimplex(face_num)
    if adj_tet is None:
        return None
    adj_face = tet.adjacentGluing(face_num)[face_num]

    return (adj_tet.index(), adj_face)

def tet_top_bottom_edge_nums(tet_vert_coorientations, tet_num):
    top_vert_nums = get_tet_top_vert_nums(tet_vert_coorientations, tet_num)
    bottom_vert_nums = list({0, 1, 2, 3} - set(top_vert_nums))
    top_edge_num = vert_pair_to_edge_num[tuple(sorted(top_vert_nums))]
    bottom_edge_num = vert_pair_to_edge_num[tuple(sorted(bottom_vert_nums))]
    return (top_edge_num, bottom_edge_num)

def tet_pair_at_edge(tri, edge_num, coors, below=False):
    edge = tri.edge(edge_num)
    for embed in edge.embeddings():
        tet = embed.tetrahedron()
        top_edge_num, bottom_edge_num = tet_top_bottom_edge_nums(coors, tet.index())
        if below:
            if tet.edge(top_edge_num) == edge:
                return (tet.index(), top_edge_num)
        elif tet.edge(bottom_edge_num) == edge:
            return (tet.index(), bottom_edge_num)
        
    return None

def face_num_meeting_edge(tri, tet_num, edge_num, coors, top=True):
    tet = tri.simplex(tet_num)
    mapping = tet.edgeMapping(edge_num)
    meeting = [f for f in range(4) if f not in (mapping[0], mapping[1])]
    top_verts = set(get_tet_top_vert_nums(coors, tet_num))
    allowed = set(range(4)) - top_verts if top else top_verts
    candidates = [f for f in meeting if f in allowed]
    assert len(candidates) == 1  # exactly one face meets this edge

    return candidates[0]

def get_top_pair(face, coors):
    for embed in face.embeddings():
        tet_num = embed.simplex().index()
        face_num = embed.vertices()[3]
        if coors[tet_num][face_num] == -1:
            return (tet_num, face_num)
    return None

def edge_sides_w_diags(tri, angle, coors):
    face_circles, tet_circles = edge_side_circle(tri, return_tets=True)

    out = []
    for edge_num in range(tri.countEdges()):
        fan = tet_circles[edge_num]
        pi_positions = []
        for i, tet_side in enumerate(fan):
            tet_num, edge_in_tet = tet_side
            if edge_in_tet in (angle[tet_num], 5 - angle[tet_num]):
                pi_positions.append(i)
        assert len(pi_positions) == 2  # taut edge has two pi sides
        p0, p1 = pi_positions
        pi_a, pi_b = (fan[p0], fan[p1])
        side_A = [fan[i] for i in range(p0 + 1, p1)]
        side_B = [fan[i] for i in range(p1 + 1, len(fan))] + [fan[i] for i in range(0, p0)]
        below_pair = tet_pair_at_edge(tri, edge_num, coors, below=True)
        if pi_a == below_pair:
            pi_below, pi_above = (pi_a, pi_b)
        elif pi_b == below_pair:
            pi_below, pi_above = (pi_b, pi_a)
        else:
            pi_below, pi_above = (pi_a, pi_b)

        out.append((side_A, [pi_below], side_B, [pi_above]))
        
    return out

def extreme_tet_side_in_ball(ball, edge_num, parity='up'):
    proj_num = ball.ledge_num_to_proj_num(edge_num)
    if parity == 'up':
        extreme = tet_pair_at_edge(ball.base_tri, proj_num, ball.base_coors)
    else:
        extreme = tet_pair_at_edge(ball.base_tri, proj_num, ball.base_coors, below=True)

    mapping = ball.edge_sides_in_quotient_mapping(edge_num)
    if extreme in mapping:
        return mapping[extreme]
    
    embed = ball.tri.edge(edge_num).embedding(0)

    return embedge_to_tet_num_pair(embed)

def equatorial_boundary_edges(ball, parity='up'):
    # The set of ball edges on the boundary of the upper, parity='up', or lower, parity='down', disk 
    
    tri = ball.tri
    disk_face_indices = set(ball.get_upper_faces(parity=parity))
    edge_face_count = {}
    for face_index in disk_face_indices:
        triangle = tri.triangle(face_index)
        for i in range(3):
            e = triangle.edge(i).index()
            edge_face_count[e] = edge_face_count.get(e, 0) + 1

    return {e for e, c in edge_face_count.items() if c == 1}

def build_equatorial_circle(ball, parity='up'):
    # Ordered (tet, edge) tet-sides tracing the ball's equator as a circle
    tri = ball.tri
    boundary_edges = equatorial_boundary_edges(ball, parity=parity)
    assert boundary_edges  # disk boundary is nonempty
    vert_to_edges = {}
    for e in boundary_edges:
        edge = tri.edge(e)
        for k in range(2):
            v = edge.vertex(k).index()
            vert_to_edges.setdefault(v, []).append(e)

    for v, es in vert_to_edges.items():
        assert len(es) == 2  # equator is a cycle

    ordered_edges = []
    start_edge = min(boundary_edges)
    prev_vert = tri.edge(start_edge).vertex(0).index()
    current_edge = start_edge

    while True:
        ordered_edges.append(current_edge)
        edge = tri.edge(current_edge)
        v0, v1 = (edge.vertex(0).index(), edge.vertex(1).index())
        next_vert = v1 if prev_vert == v0 else v0
        e0, e1 = vert_to_edges[next_vert]
        next_edge = e1 if current_edge == e0 else e0
        prev_vert, current_edge = (next_vert, next_edge)
        if current_edge == start_edge:
            break

    assert len(ordered_edges) == len(boundary_edges)  # walk covers the equator
    return [extreme_tet_side_in_ball(ball, e, parity=parity) for e in ordered_edges]

def circle_agrees_with_reference(circle, ref):
    # Does cyclic vertex list circle contain the triple ref = (a, b, c) in the same cyclic order?
    pos = {v: i for i, v in enumerate(circle)}
    a, b, c = ref
    if a not in pos or b not in pos or c not in pos:
        return True
    
    ia = pos[a]
    n = len(circle)
    db = (pos[b] - ia) % n
    dc = (pos[c] - ia) % n
    return db < dc