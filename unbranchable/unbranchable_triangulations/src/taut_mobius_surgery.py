#
# taut_mobius_surgery.py
# for performing the generalisation of veering dehn surgery on mobius band triangles
# see also veering's implementation of the same operation on veering tris

import regina
from veering.taut import is_taut, edge_num_to_vert_pair, vert_pair_to_edge_num
from veering.transverse_taut import edge_side_face_collections

def unjoin_face(tri,f):
    embed = tri.triangle(f).embedding(0)
    embed.simplex().unjoin(embed.face())
    return tri
    
def get_mobius_strip_indices(triangulation):
    return [f.index() for f in triangulation.triangles() if f.isMobiusBand()]

# here, ang = angle[tet_num], vert = num of vertex in tet 
def shares_pi_with(ang, vert):
    other_part = [1, 2, 3]
    other_part.remove(ang + 1)
    for part in ([0, ang + 1], other_part):
        if vert in part:
            part.remove(vert)
            return part.pop()

def band_vertex(face):
    # The face vertex on both copies of the identified edge of the Mobius band
    edges = [face.edge(i) for i in range(3)]
    for j in range(3):
        if edges[j] == edges[(j + 1) % 3]:
            return (j + 2) % 3
    return None

def pi_is_sharp(tri,angle,face):
    band_edges = [i for i in range(3) if i!= band_vertex(face)]
    edge_num = face.edge(band_edges[0]).index()
    pairs = edge_side_face_collections(tri,angle)[edge_num]
    fans = set()
    for i in range(2):
        pair = (face.index(),band_edges[i])
        fan = next(j for j in range(2) if pair in pairs[j])
        fans.add(fan)
    return len(fans)==1

def pi_edges_in_face(embedding, ang):
    # embedding of face, ang = angle of respective tet
    # returns which vertices are the pi edge in the face from an embedding
    verts = embedding.vertices()
    face_verts = {verts[i] for i in range(3)}
    pi_pair_1 = set(edge_num_to_vert_pair[ang])
    pi_pair_2 = set(edge_num_to_vert_pair[5-ang])

    if pi_pair_1 <= face_verts:
        in_face = pi_pair_1
    else:
        assert pi_pair_2 <= face_verts # triangulation was not taut (or even z/2z taut!): face has no pi
        in_face = pi_pair_2

    return frozenset(verts.inverse()[x] for x in in_face)


def in_face_pi_kinds(face, angle):
    # are they base or band
    B = band_vertex(face)
    if B is None:
        return None, None
    
    rest = [v for v in range(3) if v != B]
    base = frozenset(rest)
    kinds = []
    
    for k in (0, 1):
        pi_in_face = pi_edges_in_face(face.embedding(k),
                                      angle[face.embedding(k).tetrahedron().index()])
        
        if pi_in_face == base:
            kinds.append("base")

        else:
            kinds.append(pi_in_face)

    return tuple(kinds)

def is_off_band(face, angle):
    # are they off the band
    if band_vertex(face) is None:
        return True
    k0, k1 = in_face_pi_kinds(face, angle)
    return k0 == "base" or k1 == "base" or k0 is None or k1 is None

def pis_coincide(face, angle):
    k0, k1 = in_face_pi_kinds(face, angle)
    return k0 is not None and k0 != "base" and k0 == k1

def is_surgery_admissible(tri, angle_struct, face_num):
    """
    Checks if the triangle satisfies veering pis+is sharp
    """
    # REMARK: It may be possible that you can do taut mobius surgery only using off_band + no pis_coincide...
    # That alone does not imply sharp
    face = tri.triangle(face_num)
    if not face.isMobiusBand():
        return "not_mobius"
    off = is_off_band(face, list(angle_struct))
    same = pis_coincide(face, list(angle_struct))
    sharp = pi_is_sharp(tri,angle_struct,face)
    if off or same or not sharp:
        return False
    #     return "both"
    # if off:
    #     return "off_band"
    # if same:
    #     return "coinciding_pi"
    return True

def surgery_admissible_faces(tri,angle):
    return [f for f in range(tri.countTriangles()) if is_surgery_admissible(tri, angle, f)]

def local_edge_colour(ang, edge_num):
    """
    The colours of a veering tetrahedron's edges. We use this locally about the band face
    """
    L_pair = (ang + 1) % 3
    R_pair = (ang + 2) % 3
    if edge_num in (L_pair, 5 - L_pair):
        return "blue"
    if edge_num in (R_pair, 5 - R_pair):
        return "red"
    return "pi"

def face_colour(triangulation, angle_struct, face_num):
    """
    The local veering colour of a Mobius band
    """
    tri = regina.Triangulation3(triangulation)
    angle = list(angle_struct)
    face = tri.triangle(face_num)
    if not face.isMobiusBand():
        return None

    slot_colours = []
    for x, y in ((0, 1), (0, 2), (1, 2)):
        colour = None
        for side in (0, 1):
            emb = face.embedding(side)
            tet = emb.tetrahedron()
            V = emb.vertices()
            edge_num = vert_pair_to_edge_num[tuple(sorted([V[x], V[y]]))]
            c = local_edge_colour(angle[tet.index()], edge_num)
            if c != "pi":
                colour = c
        slot_colours.append(colour)

    if None in slot_colours:
        return None
    
    reds = slot_colours.count("red")
    blues = slot_colours.count("blue")

    if reds == blues:
        return None # =>input was not surgery admissible
    elif reds>blues:
        return "red"
    else:
        return "blue"

# Now we have the setup where we can just rip veering_dehn_surgery.py and remove references to global
# veering structure. They give an explanation there. 
def taut_mobius_dehn_surgery(triangulation, angle_struct, face_num):
    angle = angle_struct[:]

    # The band must be admissible for surgery.
    if not is_surgery_admissible(triangulation, angle_struct, face_num):
        return False, False, False

    # Determine the colour
    colour = face_colour(triangulation, angle, face_num)
    assert colour is not None
    assert colour in ("red", "blue")

    tri = regina.Triangulation3(triangulation)
    face = tri.triangle(face_num)

    # Band labels from veering
    #     (a, b, c, d) on tet0
    #     (p, q, r, s) on tet1

    edges = [face.edge(i) for i in range(3)]

    B = None
    for j in range(3):
        if edges[j] == edges[(j + 1) % 3]:
            B = (j + 2) % 3
            break

    assert B is not None

    embed0 = face.embedding(0)
    embed1 = face.embedding(1)

    tet0 = embed0.tetrahedron()
    tet1 = embed1.tetrahedron()

    i0 = tet0.index()
    i1 = tet1.index()

    b = embed0.vertices()[B]
    c = shares_pi_with(angle[i0], b)
    d = embed0.vertices()[3]

    rest0 = [i for i in range(4) if i not in (b, c, d)]

    q = embed1.vertices()[B]
    r = shares_pi_with(angle[i1], q)
    s = embed1.vertices()[3]

    rest1 = [i for i in range(4) if i not in (q, r, s)]

    assert c != d
    assert r != s
    assert len(rest0) == 1 and len(rest1) == 1

    a = rest0[0]
    p = rest1[0]

    # Pop open the face
    tet0 = tri.tetrahedron(i0)
    tet0.unjoin(d)

    # Refresh tetrahedron references after changing the triangulation.
    tet0 = tri.tetrahedron(i0)
    tet1 = tri.tetrahedron(i1)

    # Now insert
    tet_new = tri.newTetrahedron()

    if colour == "red":
        plan = [
            (0, tet_new, regina.Perm4(3, 0, 1, 2)),
            (1, tet0,    regina.Perm4(c, d, a, b)),
            (2, tet1,    regina.Perm4(q, p, s, r)),
        ]

    else:  # blue
        plan = [
            (1, tet_new, regina.Perm4(1, 3, 0, 2)),
            (2, tet0,    regina.Perm4(a, b, d, c)),
            (0, tet1,    regina.Perm4(s, r, p, q)),]

    for src_facet, dst, perm in plan:
        img = perm[src_facet]

        assert dst.adjacentSimplex(img) is None

        tet_new.join(src_facet, dst, perm)

    assert tri.isValid()
    assert tri.countTetrahedra() == tri.countEdges()

    # Now determine the angle. The choice here is unique, we find it by just looking
    for new_pi in (0, 1, 2):
        candidate = angle + [new_pi]

        if is_taut(tri, candidate):
            return tri, candidate, tet_new.triangle(3).index()

    assert False

def taut_n_mobius_dehn_surgery(tri, angle, face_num, n):
    for _ in range(n):
        tri, angle, face_num = taut_mobius_dehn_surgery(tri, angle, face_num)
    return tri, angle, face_num