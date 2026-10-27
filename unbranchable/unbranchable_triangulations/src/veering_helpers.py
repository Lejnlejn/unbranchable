# veering_helpers.py

"""
We need certian functions from the veering package but in contexts with both boundary and with Z2-taut instead of Z-taut.
The altered version are here. 
"""

import regina
from copy import deepcopy

from veering.taut import apply_isom_to_angle_struct_list, edge_num_to_vert_pair, pi_edgepair, vert_pair_to_edge_num
from veering.transverse_taut import convert_tetrahedron_coorientations_to_faces, vertexSplit

def taut_regina_angle_struct_to_taut_struct(regina_angle_struct):
    """
    Convert a taut regina angle structure to list of which edge pair
    is the pi pair.
    """
    out = []
    for tet_num in range(regina_angle_struct.triangulation().countTetrahedra()):
        out.append(pi_edgepair(regina_angle_struct, tet_num))
    return out

def get_taut_structures(tri):
    """
    A wrapper for the regina function.
    """
    regina_angles = regina.AngleStructures(tri, True)
    return [taut_regina_angle_struct_to_taut_struct(angle) for angle in regina_angles]

edge_pair_to_edge_numbers_dict = {0:(0,5),1:(1,4),2:(2,3)}

def is_z2_taut(tri, angle):
    totals = [0] * tri.countEdges()
    for i, tet in enumerate(tri.tetrahedra()):
        edge_nums = edge_pair_to_edge_numbers_dict[angle[i]]
        for e in edge_nums:
            totals[tet.edge(e).index()] += 1
    # print totals
    totals = [totals[i] for i in range(tri.countEdges()) if not tri.edge(i).isBoundary()]
    return all(total%2 == 0 for total in totals) # if e is boundary, then nothing obstructs z2 tautness from e's POV


def find_z2_taut_structures(tri):
    """try all possibilities of flat tetrahedra, looking for a z2 taut structure on tri"""

    num_edges = tri.countEdges()
    remaining_edge_degrees = [e.degree() for e in tri.edges()]
    pass_back_structures = []
    try_build_z2_taut_struct(tri, [], [0]*num_edges, remaining_edge_degrees, pass_back_structures)
    return pass_back_structures

def num_z2_taut_automorphisms(tri, angle):
    """
    [LH: tinkered from taut.py]
    Finds the number of z2-taut automorphisms
    """
    all_isoms = tri.findAllIsomorphisms(tri)
    all_angles = []
    for isom in all_isoms:
        all_angles.append( apply_isom_to_angle_struct_list(angle, isom) )
    return all_angles.count(all_angles[0])


def try_build_z2_taut_struct(tri, partial_taut_structure, num_pis_at_edges, remaining_edge_degrees, pass_back_structures):
    """recursive function explores tree of partial taut structures looking for z2 taut. builds taut structure
    by adding tetrahedra in sequentially by tetrahedron number"""

    if len(partial_taut_structure) == tri.countTetrahedra():
        # print partial_taut_structure, num_pis_at_edges, veering_directions
        assert is_z2_taut(tri, partial_taut_structure)
        pass_back_structures.append(partial_taut_structure)
        return True
    #try to add next tetrahedron in
    tet_index = len(partial_taut_structure)
    tet = tri.tetrahedron(tet_index)
    for edgepair in range(3):
        try_add_z2_taut_tet(tri, partial_taut_structure, num_pis_at_edges, remaining_edge_degrees, tet, edgepair, pass_back_structures)


def try_add_z2_taut_tet(tri, partial_taut_structure, num_pis_at_edges, remaining_edge_degrees, tet, edgepair, pass_back_structures):
    """try to add the given tet, tet_index to the structure with given edgepair being pi"""
    #first see if num_pis_at_edges is ok
    new_num_pis_at_edges = deepcopy(num_pis_at_edges)
    new_remaining_edge_degrees = deepcopy(remaining_edge_degrees)

    pi_edge_1, pi_edge_2 = edgepair, 5 - edgepair
    edge_1_index = tet.edge(pi_edge_1).index()
    edge_2_index = tet.edge(pi_edge_2).index()
    new_num_pis_at_edges[edge_1_index] += 1
    new_num_pis_at_edges[edge_2_index] += 1
    for e in range(6):
        edge = tet.edge(e)
        e_index = edge.index()
        new_remaining_edge_degrees[e_index] -= 1
        assert new_remaining_edge_degrees[e_index] >= 0
        if new_remaining_edge_degrees[e_index] == 0:
            if not edge.isBoundary() and new_num_pis_at_edges[e_index] % 2 != 0:
                return False

    #survived this far, so go deeper into the tree
    new_partial_taut_structure = deepcopy(partial_taut_structure)
    new_partial_taut_structure.append(edgepair)
    try_build_z2_taut_struct(tri, new_partial_taut_structure, new_num_pis_at_edges, new_remaining_edge_degrees, pass_back_structures)
    
def reverse_tet_orientation(triangulation, tet, pi_location, merge=None):
    # only need triangulation here for diagnostics
    """
    Reglues tet into triangulation with reversed orientation.  Returns
    which permutation of the tetrahedron was used.
    """
    swaps = {0: regina.Perm4(1, 0, 2, 3), 1: regina.Perm4(2, 1, 0, 3), 2: regina.Perm4(3, 1, 2, 0)}
    swap = swaps[pi_location]
    adjtets = []
    oldadjgluings = []
    adjgluings = []
    self_faces = []
    other_faces = []
    bdry_faces=[]
    for face in range(4):
        if tet.triangle(face).isBoundary():
            bdry_faces.append(face)
            adjtets.append(None)
            adjgluings.append(None)
        else:
            adjtet = tet.adjacentTetrahedron(face)
            adjtets.append(adjtet)
            oldadjgluing = tet.adjacentGluing(face)
            oldadjgluings.append(tet.adjacentGluing(face))
            if adjtet != tet:
                adjgluings.append(oldadjgluing * swap)  # swaps 0 with 1 in this tet, not the other tet
                other_faces.append(face)
            else:
                adjgluings.append(swap * oldadjgluing * swap)
                self_faces.append(face)
    tet.isolate()  # undo all gluings
    for face in other_faces:
        tet.join(swap[face], adjtets[face], adjgluings[face])
    if len(self_faces) > 0:  # assume it must be two
        assert len(self_faces) == 2, "4 self faces, must be an error..."
        face = self_faces[0]  # only need to make one of the identifications
        tet.join(swap[face], adjtets[face], adjgluings[face])
    if merge!=None:
        # sadly, this /will/ mess up where the front edge is...
        front_edge = merge[tet.index()]
        front_verts = edge_num_to_vert_pair[front_edge]
        new_front_edge = vert_pair_to_edge_num[tuple(swap[v] for v in front_verts)]
        merge[tet.index()] = new_front_edge
    return swap

def find_orientations(triangulation):
    """
    Assuming the triangulation is orientable, finds the current
    orientations of the tetrahedra.
    """
    orientations = [0] * triangulation.countTetrahedra()
    orientations[0] = +1  # fix the orientation for tet 0
    frontier = [(0, 0), (0, 1), (0, 2), (0, 3)]
    while len(frontier) > 0:
        tet_num, face_num = frontier.pop()
        tet = triangulation.tetrahedron(tet_num)
        adjtet, adjgluing = tet.adjacentTetrahedron(face_num), tet.adjacentGluing(face_num)
        if adjtet!=None:
            if orientations[adjtet.index()] != 0:  # we have already been to that tetrahedron
                assert orientations[adjtet.index()] == adjgluing.sign() * -1 * orientations[tet.index()], "triangulation not orientable!"
                frontier.remove((adjtet.index(), adjgluing[face_num]))  # remove the reverse gluing from the frontier
            else:
                orientations[adjtet.index()] = adjgluing.sign() * -1 * orientations[tet.index()]
                for i in range(4):
                    if i != adjgluing[face_num]:  # don't go backwards
                        frontier.append((adjtet.index(), i))
    for orient in orientations:
        assert orient != 0, "triangulation not connected!"
    return orientations
    
def fix_orientations(tri, angle, merge=None, return_isom = False):
    """
    Fix the orientations of the tetrahedra in triangulation so that
    they are consistently oriented.  We choose how to flip each
    tetrahedron so that the angle structure list does not change.  If
    return_isom, return the isomorphism from the original tri to the
    fixed orientations tri
    """
    orig_tri = regina.Triangulation3(tri)
    old_orientations = find_orientations(tri)
    swaps = []
    for i, orientation in enumerate(old_orientations):
        if orientation == -1:
            swaps.append( reverse_tet_orientation(tri, tri.tetrahedron(i), angle[i], merge) )
        else:
            swaps.append( regina.Perm4() ) ## identity

    if return_isom:  
        out_isom = regina.Isomorphism3.identity(len(swaps))
        for i, p in enumerate(swaps):
            out_isom.setFacetPerm(i, p)
        return out_isom

def is_transverse_taut(triangulation, taut_angle_struct, return_type = "boolean"):
    ''' Ripped and changed from veering.transverse_taut to work for Z/2'''
    assert is_z2_taut(triangulation, taut_angle_struct)
    assert triangulation.isOriented()
    coorientations = []
    for i in range(triangulation.countTetrahedra()):
        coorientations.append([0,0,0,0]) # +1 is out of tet, -1 is into tet
    # first decide coorientation for tet 0
    edgepair = taut_angle_struct[0]
    if not add_coors_with_check(triangulation, coorientations, 0, edgepair, +1): # choose arbitrary direction for tet 0
        return False # tet 0 could have self gluings incompatible with the transverse taut structure, so we could fail here

    explored_tetrahedra = [0]
    frontier_tet_faces = [(0, 0), (0, 1), (0, 2), (0, 3)]
    while len(frontier_tet_faces) > 0:
        my_tet_num, my_face_num = frontier_tet_faces.pop()
        my_tet = triangulation.tetrahedron(my_tet_num)
        face = my_tet.triangle(my_face_num)
        if not face.isBoundary():
            neighbour_tet = my_tet.adjacentSimplex(my_face_num)
            neighbour_face_num = my_tet.adjacentGluing(my_face_num)[my_face_num]
            if neighbour_tet.index() in explored_tetrahedra:
                frontier_tet_faces.remove((neighbour_tet.index(), neighbour_face_num))
            else:
                newcoor = -coorientations[my_tet_num][my_face_num] # opposite of the adjacent face's coorientation
                edgepair = taut_angle_struct[neighbour_tet.index()]
                # now calculate the direction that the flow points through neighbour_tet
                if neighbour_face_num in [vertexSplit[edgepair][i] for i in range(2)]:
                    direction = newcoor # agrees with newcoor
                else:
                    direction = -newcoor
                if not add_coors_with_check(triangulation, coorientations, neighbour_tet.index(), taut_angle_struct[neighbour_tet.index()], direction):
                    return False
    
                frontier_tet_faces.extend( [(neighbour_tet.index(), i) for i in range(4) if i != neighbour_face_num] )
                explored_tetrahedra.append(neighbour_tet.index())
    if return_type == "tet_vert_coorientations":
        return coorientations
    elif return_type == "face_coorientations":
        return convert_tetrahedron_coorientations_to_faces(triangulation, coorientations)
    else:
        return True


def add_coors_with_check(triangulation, coorientations, tet_num, edgepair, direction):
    """
    Adds coorientation directions to coorientations for this tet_num,
    checks if they match neighbours' coorientations
    """
    # edgepair is the taut angle structure pi angle direction is +1 if
    # from the 0,1 faces to the 2,3 faces if we are in vertexSplit[0],
    # and similarly for the other cases
    for i in range(0, 2):
        if coorientations[tet_num][vertexSplit[edgepair][i]] == -direction:
            return False
        else:
            coorientations[tet_num][vertexSplit[edgepair][i]] = direction
    for i in range(2, 4):
        if coorientations[tet_num][vertexSplit[edgepair][i]] == direction:
            return False
        else:
            coorientations[tet_num][vertexSplit[edgepair][i]] = -direction
    # now check whether these new coorientations match the tetrahedra
    # these are glued to:
    tet = triangulation.tetrahedron(tet_num)
    for face in range(4):
        if not tet.triangle(face).isBoundary():
            adjtet, adjgluing = tet.adjacentTetrahedron(face), tet.adjacentGluing(face)
            adjtet_num = adjtet.index()
            if coorientations[tet_num][face] * coorientations[adjtet_num][adjgluing[face]] == 1:
                return False # failed neighbours check
    return True # passed neighbours checks

def is_taut(tri, angle, return_totals = False):
    ''' adapted for the case where the edge is boundary (!) '''
    #n = tri.countTetrahedra()
    #assert n == tri.countEdges()
    #assert len(angle) == n
    non_bdry_edges,bdry_edges=[],[]
    n =  tri.countEdges()
    for i in range(n):
        if tri.edge(i).isBoundary():
            bdry_edges.append(i)
        else:
            non_bdry_edges.append(i)
    totals = [0] * n
    for i, tet in enumerate(tri.tetrahedra()):
        edge_nums = edge_pair_to_edge_numbers_dict[angle[i]]
        for e in edge_nums:
            totals[tet.edge(e).index()] += 1
    if return_totals == True:
        return totals
    
    return all(totals[i] == 2 for i in non_bdry_edges) and all(totals[i] <= 2 for i in bdry_edges)

def isosig_to_tri_angle(isosig, return_isom = False):
    """
    Modified from veering.taut to allow preangle structures. That is, any local assignment of taut angles to each tetrahedron.
    Given a taut isosig, returns an oriented regina triangulation and
    the list of angles for the taut angle structure, for the new
    labelling.  If return_isom, returns the isomorphism from the
    regina isosig triang to the oriented triangulation.
    """
    data = isosig.split("_")
    isosig, angle = data[0], data[1]
    tri = regina.Triangulation3.fromIsoSig(isosig)
    angle = [int(i) for i in list(angle)]
    isom = fix_orientations(tri, angle, return_isom = return_isom)  # this does not alter what the angles_string should be
    assert tri.isOriented()
    #assert is_taut(tri, angle) # disable for generalised taut structures
    if return_isom:
        return tri, angle, isom
    else:
        return tri, angle
