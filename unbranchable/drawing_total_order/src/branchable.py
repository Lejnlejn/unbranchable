# branchable.py

# This is a local copy to keep the remaining features of this project isolated.
# For working with branchings on triangulations. This breaks down to "Z/2Z taut merging structures",
# which are `the' natural dual objects to branchings. t
# A branching structure is encoded by two tuples: a list of pi pair for each tetrahedron, and a list of which edge in the tetrahedron
# is a front edge. This specifies a unique smoothing on the dual spine of the tetrahedron.

import regina
from veering.taut import apply_isom_to_angle_struct_list, edge_num_to_vert_pair, vert_pair_to_edge_num, vert_pair_to_edge_pair #, edge_pair_to_edge_numbers #,find_orientations #,taut_regina_angle_struct_to_taut_struct, get_taut_structures
from veering.transverse_taut import edge_side_face_collections, convert_tetrahedron_coorientations_to_faces
# from veering.z2_taut import find_z2_taut_structures #, is_z2_taut is adapted below to the with-boundary case. find_z2 should also be.
from copy import deepcopy
from veering.taut import pi_edgepair, isosig_from_tri_angle
from veering.transverse_taut import vertexSplit


# Small Regina-convention helpers used throughout this module.
edge_pair_to_vert_pairs = {
    0: ((0, 1), (2, 3)),
    1: ((0, 2), (1, 3)),
    2: ((0, 3), (1, 2)),
}


def get_fembed(tri, pair):
    tet_num, face_num = pair
    return next(
        f
        for f in tri.simplex(tet_num).triangle(face_num).embeddings()
        if f.face() == face_num and f.simplex().index() == tet_num
    )

# Unfortunately, to any dear reader, we currently need to slightly adjust various functions from veering
# for our purpose. So, we do this first.

# Veering rips
#================================================================================================
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
#================================================================================================
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


# Back to branchings
#================================================================================================
def load_tri_angle_merge_from_sig(isosig, return_isom = False, apply_isom=None):
    """
    LH: Modified from veering.taut to allow preangle structures.
    Given a taut isosig, returns an oriented regina triangulation and
    the list of angles for the taut angle structure, for the new
    labelling. 
    
    If return_isom, returns the isomorphism from the
    regina isosig triang to the oriented triangulation.

    If apply_isom = regina.isom3, applies isom to tri,angle,merge

    If the sig is only a taut isoSig, it will try to find a transitive merging structure, and otherwise will return an error. 
    """
    data = isosig.split("_")
    if len(data)==3:
        isosig, angle, merge = data
        tri = regina.Triangulation3.fromIsoSig(isosig)
        angle = [int(i) for i in list(angle)]
    else:
        isosig, angle = data
        tri = regina.Triangulation3.fromIsoSig(isosig)
        angle = [int(i) for i in list(angle)]
        merges = find_transitive_z2_merge_structures_fast(tri,angle) 
        if len(merges)==1:
            merge = merges[0]
        else:
            raise ValueError("Taut structure does not have a transitive merging structure")
    
    merge = [int(m) for m in list(merge)]
    isom = fix_orientations(tri, angle, merge, return_isom = return_isom)  # this does not alter what the angles_string should be
    assert tri.isOriented()
    #assert is_taut(tri, angle) # disable for generalised taut structures
    if apply_isom!=None:
        og_isom = apply_isom.inverse()
        tri = og_isom(tri)
        apply_isom_to_angle_struct_list(angle,og_isom)
        apply_isom_to_front_edges(merge,og_isom)
  
    if return_isom:
        return tri, angle, merge, isom
    else:
        return tri, angle, merge

#========
def load_both_coor_choices(edge_split):
    ''' given an `edge splitting' loads the two alternating assignments of 0,1s '''
    coor_choices=[]
    if len(edge_split)==1:
        coor_choices=[[0],[1]]
    else:
        edge_choices=[]
        assert len(edge_split)%2==0 #should be even, else not Z2 taut
        for j in range(2):
            choice = [j for k in edge_split]
            swap = (j+1)%2
            for i in range(int(len(edge_split)/2)):
                choice[2*i]=swap
            edge_choices.append(choice)
        coor_choices=edge_choices
    return coor_choices


def face_num_to_edges(face_num):
    face_verts = list(set(range(4))-{face_num})
    vert_pairs = [(face_verts[0],face_verts[1]),(face_verts[0],face_verts[2]),(face_verts[1],face_verts[2])]
    return list(vert_pair_to_edge_num[pair] for pair in vert_pairs)

# load this up so its only ever a call-up again
face_num_to_edges_dict = {i: face_num_to_edges(i) for i in range(4)}

def face_num_to_edges(face_num):
    return face_num_to_edges_dict[face_num]

def face_num_to_large_edge(face,tri,angle,merge):
    fembed = tri.triangle(face).embedding(0) # assumes consistent!
    tet_num =fembed.simplex().index()
    internal_face = fembed.face()
    edges = face_num_to_edges_dict[internal_face]
    front_edge = merge[tet_num]
    if front_edge in edges:
        large_edge = front_edge
    else:
        # its whatever your pi edge is
        pi_edges = edge_pair_to_edge_numbers_dict[angle[tet_num]]
        large_edge = next(p for p in pi_edges if p in edges)
    
    hit_verts = [fembed.vertices().inverse()[i] for i in edge_num_to_vert_pair[large_edge]]
    missing_vert = next(i for i in range(3) if i not in hit_verts) # this is the image edge
    return missing_vert
    
def merge_to_face_large_edges(tri,angle,merge):
    output=[]
    for face in range(tri.countTriangles()):
        output.append(face_num_to_large_edge(face,tri,angle,merge))
    return output
    
def merge_to_coors(tri,angle,merge):
    '''assumes that merge is indeed a valid merge structure'''
    face_large_edges = merge_to_face_large_edges(tri,angle,merge)
    edge_split = edge_sides_splitted(tri,angle)
    # now for each edge_split, we will want to check whether or not that is a large edge or not.
    output=[]
    for edge in edge_split:
        edge_put=[]
        for split in edge:
            face_num,edge_num = split[0] # they will all be the same
            coor = int(face_large_edges[face_num]==edge_num) # i.e. 1 = large edge is how we meet the edge = 'inward'
            edge_put.append(coor)
        output.append(edge_put)
    return output

def coor_to_mapping(coor,tri,angle):
    edge_side_chunks = edge_sides_splitted(tri,angle)
    mapping={}
    for e in range(tri.countEdges()):
        edge_coors = coor[e]
        chunks = edge_side_chunks[e]
        for i in range(len(chunks)):
            local_coor = edge_coors[i]
            chunk = chunks[i]
            for c in chunk:
                mapping.update({c:local_coor})
    return mapping
    
def is_merge_mapping(mapping,tri):
    '''Determines if a mapping (face,edge): in or out determines a merge structure'''
    for t in range(tri.countTriangles()):
        tri_sides = [(t,i) for i in range(3)]
        tally = 0
        for pair in tri_sides:
            tally = tally+mapping[pair]
        if tally==1: # this forces the convention that into face=0 and out of face=1
            pass
        else:
            return False
    return True

def is_z2_merge_coor(coor,tri,angle):
    return is_merge_mapping(coor_to_mapping(coor,tri,angle),tri)
    
def isosig_to_tri_angle(isosig, return_isom = False):
    """
    LH: Modified from veering.taut to allow preangle structures.
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
    
def in_face_to_corresponding_coor(tri,edge_splits,in_pair):
    face_num,internal = in_pair
    global_edge_num = tri.triangle(face_num).edge(internal).index()
    edge_split = edge_splits[global_edge_num]
    slot = [x for x in range(len(edge_split)) if in_pair in edge_split[x]][0]
    both_coors = load_both_coor_choices(edge_split)
    matching_coor = [i for i in range(2) if both_coors[i][slot]==0][0] # the coor where this face points into the triangle
    return (global_edge_num,both_coors[matching_coor])

def update_coor_from_in_face(tri,alive_coors,edge_splits,in_pair):
    ''' returns True, updated coor if possible, otherwise, the coor is inconsistent, and returns False, None'''
    edge,coor = in_face_to_corresponding_coor(tri,edge_splits,in_pair)
    #print(in_pair)
    if edge in alive_coors.keys():
        if alive_coors[edge]==coor:
            return True, alive_coors, None
        else:
            return False, None, None
    else:
        alive_coors.update({edge:coor})
        #print(edge)
        return True, alive_coors, edge

def step_new_coors_out_one_from_a_single_coor(tri,new_edges,current_coor,edge_splits):
    ''' propagate a coor assignment to neighbouring faces'''
    next_edges=[]
    for new_edge in new_edges:
        local_coor = current_coor[new_edge]
        out_indices = [i for i in range(len(local_coor)) if local_coor[i]==1]
        out_pairs = []
        edge_split = edge_splits[new_edge]
        for i in out_indices:
            out_pairs = out_pairs+edge_split[i] 
            
        # now we have a collection of out_pairs, and we want to try and update them all
        
        for pair in out_pairs:
            face_num,internal_edge = pair
            in_pairs = [(face_num,i) for i in range(3) if i!=internal_edge]
            for in_pair in in_pairs:
                #print
                verdict, current_coor, new_edge = update_coor_from_in_face(tri,current_coor,edge_splits,in_pair)
                
                if not verdict:
                    return False, None, None
                elif new_edge==None:
                    # we didn't change anything, it was just consistent
                    pass
                else:
                    next_edges.append(new_edge)
    return True, current_coor, next_edges

def exhaust_new_coor(tri,current_coors,new_edge_num,edge_splits):
    ''' After a new coor has been added, propagate it through the triangulation until either we get an inconsistency
    or run out of new places to propagate.'''
    new_edges=list(current_coors.keys())
    next_edges = new_edges[:]
    while len(next_edges)>0:
        verdict, current_coor, next_edges = step_new_coors_out_one_from_a_single_coor(tri,new_edges,current_coors,edge_splits)
        if not verdict:
            return False 
    return current_coors

def find_z2_merge_coors(tri,angle):
    ''' Search for z2 merge coors. There are two choices of merge coors at each edge. We pick one, then propagate the consequences.
    We choose again if we have to, until all edges are assigned. '''
    edge_count=tri.countEdges()
    edge_splits = edge_sides_splitted(tri,angle)
    edge_num=0
    starting_num=0
    edge_split = edge_splits[starting_num]
    edge_coors = load_both_coor_choices(edge_split)
    alive_coors = [{starting_num:x} for x in edge_coors]
    unfinished_coors = [(x,edge_num) for x in alive_coors]
    finished_coors=[]
    while len(unfinished_coors)>0:
        coor,edge_num = unfinished_coors.pop()
        output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
                
        exhausted_coors=[]
        if output!=False:
            coor = output # the one with edges proliferated
            
            if len(coor.keys())==edge_count:
                # then we've spread everything we can, so we check if it the spreadings were consistent with being merge
                if is_z2_merge_coor(coor,tri,angle):
                    finished_coors.append(coor)
            else:
                # then coor is unfinished. We append both choices for any unassigned edge and continue
                unassigned_edges = set(range(edge_count))-set(coor.keys())
                edge_num = unassigned_edges.pop()
                edge_split = edge_splits[edge_num]
                new_coors = load_both_coor_choices(edge_split)
                for new_coor in new_coors:
                    joined_coor = deepcopy(coor)
                    joined_coor.update({edge_num:new_coor})
                    # now add the next batch of both of them in 
                    unfinished_coors.append((joined_coor,edge_num))
                    
    # now convert coors to front edges
    output=[]
    for coor in finished_coors:
        output.append(coor)
    return output
    
def has_z2_merge_coors(tri,angle):
    ''' Same as find_z2_merge_coors, just cutting short for a yes/no answer '''
    edge_count=tri.countEdges()
    edge_splits = edge_sides_splitted(tri,angle)
    edge_num=0
    starting_num=0
    edge_split = edge_splits[starting_num]
    edge_coors = load_both_coor_choices(edge_split)
    alive_coors = [{starting_num:x} for x in edge_coors]
    unfinished_coors = [(x,edge_num) for x in alive_coors]
    finished_coors=[]
    
    while len(unfinished_coors)>0:
        coor,edge_num = unfinished_coors.pop()
        output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
                
        exhausted_coors=[]
        if output!=False:
            coor = output # the one with edges proliferated, sure
            
            if len(coor.keys())==edge_count:
                # then we've spread everything, so we can just check and try
                if is_z2_merge_coor(coor,tri,angle):
                    return True
            else:
                # then coor is unfinished. We want to add a new edge to get it popping
                unassigned_edges = set(range(edge_count))-set(coor.keys())
                edge_num = unassigned_edges.pop()
                edge_split = edge_splits[edge_num]
                new_coors = load_both_coor_choices(edge_split)
                for new_coor in new_coors:
                    joined_coor = deepcopy(coor)
                    joined_coor.update({edge_num:new_coor})
                    # now add the next batch of both of them in 
                    unfinished_coors.append((joined_coor,edge_num))
                
    return False

def find_z2_merge_structures(tri,angle=None):
    ''' Returns all front edge data that form a z2-merge structure with tri, angle. If no angle is specified, returns pairs (angle,merge)
    that are z2 merge structures on tri.'''
    edge_count=tri.countEdges()
    base_angle = angle
    if angle==None:
        angles = find_z2_taut_structures(tri)
    else:
        angles = [angle]
    return_output=[]
    for angle in angles:
        edge_splits = edge_sides_splitted(tri,angle)
        edge_num=0
        starting_num=0
        edge_split = edge_splits[starting_num]
        edge_coors = load_both_coor_choices(edge_split)
        alive_coors = [{starting_num:x} for x in edge_coors]
        unfinished_coors = [(x,edge_num) for x in alive_coors]
        finished_coors=[]
        
        while len(unfinished_coors)>0:
            coor,edge_num = unfinished_coors.pop()
            output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
                    
            exhausted_coors=[]
            if output!=False:
                coor = output # the one with edges proliferated, sure
                
                if len(coor.keys())==edge_count:
                    # then we've spread everything, so we can just check and try
                    if is_z2_merge_coor(coor,tri,angle):
                        finished_coors.append(coor)
                else:
                    # then coor is unfinished. We want to add a new edge to get it popping
                    unassigned_edges = set(range(edge_count))-set(coor.keys())
                    edge_num = unassigned_edges.pop()
                    edge_split = edge_splits[edge_num]
                    new_coors = load_both_coor_choices(edge_split)
                    for new_coor in new_coors:
                        joined_coor = deepcopy(coor)
                        joined_coor.update({edge_num:new_coor})
                        # now add the next batch of both of them in 
                        unfinished_coors.append((joined_coor,edge_num))
    
        for coors in finished_coors:
            merge = coor_mapping_to_front_edges(tri,angle,coors)
            if base_angle!=None:
                return_output.append(merge)
            else:
                return_output.append((angle,merge))
        
    return return_output

def split_pairs(tri,preangle):
    ''' Given a preangle structure on tri, determine the (face,edge) pairs that bound pis'''
    output=[[] for i in range(tri.countEdges())]
    for tet_num in range(tri.countTetrahedra()):
        ang = preangle[tet_num]
        tet = tri.tetrahedron(tet_num)
        vert_pairs= edge_pair_to_vert_pairs[ang]
        
        # for each vert pair, we want to find the two faces either side of it
        #tetput=[]
        for p in vert_pairs:
            
            face_side_pair=[]
            edge_num = vert_pair_to_edge_num[p]
            global_edge_num = tet.edge(edge_num).index()
            face_nums = list(set(range(4))-set(p))
            for f in face_nums:
                fembed = get_fembed(tri,(tet_num,f))
                triangle= tet.triangle(f).index()
                verts=fembed.vertices().inverse()
                image_verts = [verts[x] for x in p]
                internal_edge = (set(range(3))-set(image_verts)).pop()
                pair = (triangle,internal_edge)
                face_side_pair.append(pair)
            output[global_edge_num].append(face_side_pair)
    
    return output

def edge_side_circle(tri,return_tets=False,return_tet_pairs=False):
    """
    Return (face,edge) pairs about an edge following in a circle around the edge
    essentialy a pare down of `edge_side_face_collections' from veering.transverse_taut
    If return tets: also returns (tet,edge) pairs around the circle
    If return_tet_pairs: also returns (tet,face) pairs around the circle
    """
    out_triangles = []
    out_tets = []
    out_tet_pairs=[]
    for edge_num in range(tri.countEdges()):
        edge = tri.edge(edge_num)
        embeddings = edge.embeddings()
        triangle_num_sets = []
        tet_num_set = []
        tet_face_set=[]
        which_set_to_add_to = 0
        for embed in embeddings:
            tet = embed.tetrahedron()
            vert_perm = embed.vertices()
            trailing_vert_num, leading_vert_num = vert_perm[2], vert_perm[3]
            f = tet.triangle(leading_vert_num)
            fmapping = tet.faceMapping(2, leading_vert_num)
            index_of_opposite_vert_in_f = fmapping.inverse()[trailing_vert_num]
            assert f.edge(index_of_opposite_vert_in_f) == edge
            triangle_num_sets.append((f.index(), index_of_opposite_vert_in_f))
            # add in the face tet data, too.
            tet_num_set.append((tet.index(), embed.face()))
            # for the pair recording
            tet_face_set.append(((tet.index(),leading_vert_num),(tet.index(),trailing_vert_num)))
            #assert 
        
        out_triangles.append(triangle_num_sets)
        out_tets.append(tet_num_set)
        tet_pair_set=[]
        for i in range(len(tet_face_set)):
            j=(i+1)%len(tet_face_set)
            pair0 = tet_face_set[i][1]
            pair1 = tet_face_set[j][0]
            if not edge.isBoundary():
                assert tri.simplex(pair0[0]).triangle(pair0[1]).index()==tri.simplex(pair1[0]).triangle(pair1[1]).index()==triangle_num_sets[j][0]
            # rotate this so that actually it is face 0 that locks in with the tet_pair
            tet_pair_set.append((pair0,pair1))
        #tet_pair_set = [tet_pair_set[-1]]+tet_pair_set[1:]
        out_tet_pairs.append(tet_pair_set)
    if not return_tets:
        return out_triangles
    elif not return_tet_pairs:
        return out_triangles, out_tets
    else:
        return out_triangles, out_tets, out_tet_pairs
        
def edge_sides_splitted(tri,angle):
    ''' Returns the edge_side_circle(tri) at each edge partitioned into sets of (face,pairs) that are not separated by any pis'''
    edge_sides = edge_side_circle(tri)
    #edge_sides = get_edge_sides_as_circle(tri,base_angle)
    splits=split_pairs(tri,angle)
    edge_chunks=[]
    for edge_num in range(tri.countEdges()):
        all_chunks=[]
        edge_side = edge_sides[edge_num]
        candidates = list(edge_side)
        current_chunk=[]
        all_chunks=[]
        splitters = splits[edge_num]
        if len(candidates)==1:
            # this should only happen for a boundary edge. And we're still getting screwed over here :(
            all_chunks = [[candidates[0]]]
        else:
            for i in range(len(candidates)):
                j = (i+1)%len(candidates)
                p0,p1 = candidates[i],candidates[j]
                current_chunk.append(p0)
                if [p0,p1] in splitters or [p1,p0] in splitters:
                    all_chunks.append(current_chunk)
                    current_chunk=[]
        
            # now we just need to update if there was NOT a split across the back
            if [p0,p1] not in splitters and [p1,p0] not in splitters and len(splitters)>0:
                #print('hi')
                if len(all_chunks)>0:
                    updated_chunk = updated_chunk = current_chunk+all_chunks[0] # careful here: used to be all_chunks[0]+current_chunk
                    all_chunks = all_chunks[1:]
                    all_chunks.append(updated_chunk)
            
        if len(splitters)==0:
            all_chunks=[edge_side]
        edge_chunks.append(all_chunks)
    return edge_chunks

def has_z2_merge_structure(tri,check_transverse=False):
    ''' Do the same as finding structures above, but just go quicker if you want a yes or no'''
    # initiliase the search
    angles = find_z2_taut_structures(tri)
    if check_transverse:
        angles = [a for a in angles if is_transverse_taut(tri,a)]
         
    for angle in angles:
        edge_count=tri.countEdges()
        edge_splits = edge_sides_splitted(tri,angle)
        edge_num=0
        edge_split = edge_splits[edge_num]
        edge_coors = load_both_coor_choices(edge_split)
        alive_coors = [{edge_num:x} for x in edge_coors]
        unfinished_coors = [(x,edge_num) for x in alive_coors]
        finished_coors=[]
        
        while len(unfinished_coors)>0:
            coor,edge_num = unfinished_coors.pop()
            output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
                    
            exhausted_coors=[]
            if output!=False:
                coor = output # the one with edges proliferated
                
                if len(coor.keys())==edge_count:
                    if is_z2_merge_coor(coor,tri,angle):
                        return True
                    else:
                        pass
                else:
                    # then coor is unfinished. We want to add a new edge to get it popping
                    unassigned_edges = set(range(edge_count))-set(coor.keys())
                    edge_num = unassigned_edges.pop()
                    edge_split = edge_splits[edge_num]
                    new_coors = load_both_coor_choices(edge_split)
                    for new_coor in new_coors:
                        joined_coor = deepcopy(coor)
                        joined_coor.update({edge_num:new_coor})
                        # now add the next batch of both of them in 
                        unfinished_coors.append((joined_coor,edge_num))
    return False
    
def has_z_merge_structure(tri,check_transverse=False,return_cert = False):
    ''' Restrict to taut angles and search for merge structures'''
    # initiliase the search
    
    angles = get_taut_structures(tri)
    if check_transverse:
        angles = [a for a in angles if is_transverse_taut(tri,a)]
    for angle in angles:
        edge_count=tri.countEdges()
        edge_splits = edge_sides_splitted(tri,angle)
        edge_num=0
        edge_split = edge_splits[edge_num]
        edge_coors = load_both_coor_choices(edge_split)
        alive_coors = [{edge_num:x} for x in edge_coors]
        unfinished_coors = [(x,edge_num) for x in alive_coors]
        finished_coors=[]
        
        while len(unfinished_coors)>0:
            coor,edge_num = unfinished_coors.pop()
            output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
                    
            exhausted_coors=[]
            if output!=False:
                coor = output # the one with edges proliferated, sure
                
                if len(coor.keys())==edge_count:
                    if is_z2_merge_coor(coor,tri,angle):
                        if return_cert:
                            return isosig_from_tri_angle(tri,angle)
                        else:
                            return True
                    else:
                        pass
                else:
                    # then coor is unfinished. We want to add a new edge to get it popping
                    unassigned_edges = set(range(edge_count))-set(coor.keys())
                    edge_num = unassigned_edges.pop()
                    edge_split = edge_splits[edge_num]
                    new_coors = load_both_coor_choices(edge_split)
                    for new_coor in new_coors:
                        joined_coor = deepcopy(coor)
                        joined_coor.update({edge_num:new_coor})
                        # now add the next batch of both of them in 
                        unfinished_coors.append((joined_coor,edge_num))
    return False

def find_z_merge_structures(tri):
    ''' Restrict to taut setting'''
    # initiliase the search
    angle_coor_pairs=[]
    angles = get_taut_structures(tri)
    final_output=[] # nah there's no way...
    for angle in angles:
        edge_count=tri.countEdges()
        edge_splits = edge_sides_splitted(tri,angle)
        edge_num=0
        edge_split = edge_splits[edge_num]
        edge_coors = load_both_coor_choices(edge_split)
        alive_coors = [{edge_num:x} for x in edge_coors]
        unfinished_coors = [(x,edge_num) for x in alive_coors]
        finished_coors=[]
        
        while len(unfinished_coors)>0:
            coor,edge_num = unfinished_coors.pop()
            output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
                    
            exhausted_coors=[]
            if output!=False:
                coor = output # the one with edges proliferated, sure
                
                if len(coor.keys())==edge_count:
                    if is_z2_merge_coor(coor,tri,angle):
                        angle_coor_pairs.append((angle,coor))
                    else:
                        pass
                else:
                    # then coor is unfinished. We want to add a new edge to get it popping
                    unassigned_edges = set(range(edge_count))-set(coor.keys())
                    edge_num = unassigned_edges.pop()
                    edge_split = edge_splits[edge_num]
                    new_coors = load_both_coor_choices(edge_split)
                    for new_coor in new_coors:
                        joined_coor = deepcopy(coor)
                        joined_coor.update({edge_num:new_coor})
                        # now add the next batch of both of them in 
                        unfinished_coors.append((joined_coor,edge_num))

        for pair in angle_coor_pairs:
            angle,coor = pair
            final_output.append((angle,coor_mapping_to_front_edges(tri,angle,coor)))
    
    return final_output

def find_merge_structures(tri,taut_over='Z'):
    ''' Enumerate all structures of the given taut type'''
    # initiliase the search
    angle_coor_pairs=[]
    return_output=[]

    if taut_over=='Z':
        angles = get_taut_structures(tri)
    else:
        assert taut_over=='Z2'
        angles = find_z2_taut_structures(tri)
    
    for angle in angles:
        edge_count=tri.countEdges()
        edge_splits = edge_sides_splitted(tri,angle)
        edge_num=0
        edge_split = edge_splits[edge_num]
        edge_coors = load_both_coor_choices(edge_split)
        alive_coors = [{edge_num:x} for x in edge_coors]
        unfinished_coors = [(x,edge_num) for x in alive_coors]
        finished_coors=[]
        
        while len(unfinished_coors)>0:
            coor,edge_num = unfinished_coors.pop()
            output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
            exhausted_coors=[]
            if output!=False:
                coor = output # the one with edges proliferated, sure
                
                if len(coor.keys())==edge_count:
                    if is_z2_merge_coor(coor,tri,angle):
                        angle_coor_pairs.append((angle,coor))
                    else:
                        pass
                else:
                    # then coor is unfinished. We want to add a new edge to get it popping
                    unassigned_edges = set(range(edge_count))-set(coor.keys())
                    edge_num = unassigned_edges.pop()
                    edge_split = edge_splits[edge_num]
                    new_coors = load_both_coor_choices(edge_split)
                    for new_coor in new_coors:
                        joined_coor = deepcopy(coor)
                        joined_coor.update({edge_num:new_coor})
                        # now add the next batch of both of them in 
                        unfinished_coors.append((joined_coor,edge_num))
        for pair in angle_coor_pairs:
            angle,coor = pair
            return_output.append((angle,coor_mapping_to_front_edges(tri,angle,coor)))
    return return_output

# with 1 being pointing into the edge, what we want to do now is look for the face pairs between the face
def face_pair_to_tet_in_between(tri,pair):
    '''Given a pair of common faces about an edge, convert it to the tet pair (tet,edge) between them
    this assumes that assume none of the edges are flat and have one tet around them, so that
    the tet is uniquely identified as the tet with embeddings in both faces
    for example, front edges->face large edges-> coors'''
    faces = [tri.triangle(x[0]) for x in pair]
    indices = [x[1] for x in pair]
    fembeds = [list(f.embeddings()) for f in faces]
    flattened = fembeds[0]+fembeds[1]
    tets_seen = [fembed.simplex().index() for fembed in flattened]
    if len(set(tets_seen))!=2:
        # then there is a shortcut
        tet_num_between = max(set(tets_seen), key=tets_seen.count)
        tet = tri.simplex(tet_num_between)
        # now we need to find the edge we're after
        edge_matches=[]
        for i in range(2):
            face_num,edge_num = pair[i]
            verts = list(set(range(3))-{edge_num})
            for fembed in fembeds[i]:
                if fembed.simplex()==tet:
                    image_verts = tuple(fembed.vertices()[x] for x in verts)
                    image_edge = vert_pair_to_edge_num[image_verts]
                    edge_matches.append(image_edge)
        edge_hit = max(set(edge_matches), key=edge_matches.count)
        return (tet_num_between,edge_hit)
    else:
        tet_nums = set(tets_seen)
        all_pairs = []
        for i in range(2):
            face_num,edge_num = pair[i]
            verts = list(set(range(3))-{edge_num})
            for fembed in fembeds[i]:
                tet_num = fembed.simplex().index()
                image_verts = tuple(fembed.vertices()[x] for x in verts)
                image_edge = vert_pair_to_edge_num[image_verts]
                all_pairs.append((tet_num,image_edge))
        return max(set(all_pairs), key=all_pairs.count)
    
def coor_mapping_to_front_edges(tri,angle,coor_mapping):
    # takes a coor mapping and returns a list whose ith index is the induced front edge of the ith tet 
    fedges=[None]*tri.countTetrahedra()
    for edge_num in range(tri.countEdges()):
        edge_side_split = edge_sides_splitted(tri,angle)[edge_num]
        #print(edge_side_split)
        parities = coor_mapping[edge_num]
        in_slots = [i for i in range(len(parities)) if parities[i]==1]
        
        for slot in in_slots:
            in_faces = edge_side_split[slot]
                
            for j in range(len(in_faces)-1):
                front_pair = (in_faces[j],in_faces[j+1])
                tet, fedge = face_pair_to_tet_in_between(tri,front_pair)
                fedges[tet] = fedge
                
            if len(edge_side_split)==1:
                # We wrap around when the edge is dual to a source disk. This does not work when the edge is flat;
                # we amend this manually
                e = tri.edge(edge_num)
                if e.degree()==2:
                    for embedge in e.embeddings():
                        tet,num = embedge.simplex().index(),embedge.face()
                        fedges[tet] = num
                else:
                    front_pair = (in_faces[-1],in_faces[0])
                    tet, fedge = face_pair_to_tet_in_between(tri,front_pair)
                    fedges[tet] = fedge
                    #print('hi',edge_num,tet)
    return fedges

def inverted_merge_angle(tri,angle,front_edges):
    return [vert_pair_to_edge_pair[edge_num_to_vert_pair[fedge]] for fedge in front_edges]

def lazy_complete_front_edges(tri,angle,partial_front_edges):
    # so without even reducing the search, you can do a cross-comparison
    # partial_front_edges has None where it does not know what the front edges are
    output=[]
    coors = find_z2_merge_coors(tri,angle)
    for coor in coors:
        front_edges = coor_mapping_to_front_edges(tri,angle,coor)
        check=True
        for i in range(len(partial_front_edges)):
            if partial_front_edges[i]==front_edges[i] or partial_front_edges[i]==None:
                pass
            else:
                check=False
                break
        if check:
            output.append(front_edges)

    return output
    
def is_z_merge(tri,angle,front_edges):
    ''' Determines if tri, angle, front edges determine a taut merging structure '''
    # tri should be taut
    if not is_taut(tri,angle):
        return False

    # front_edges should not be where the angles are
    for i in range(len(angle)):
        if front_edges[i] in {angle[i],5-angle[i]}:
            return False

    # now check if everything matches over faces
    # every face-pairing preserved the front edge assigned locally by pis and front-edges 
    for face in tri.triangles():
        front_edges_in_internal_face=set()
        for fembed  in face.embeddings():
            face_num = fembed.face()
            tet_num = fembed.simplex().index()
            tet_fedge = front_edges[tet_num]
            ang = angle[tet_num]
            pis = [ang,5-ang]
            # need to determine if it is a pi or the front edge that endows the front, here
            face_edges = face_num_to_edges_dict[face_num]
            local_front=None
            if tet_fedge in face_edges:
                local_front = tet_fedge
            else:
                local_front = set(pis).intersection(set(face_edges)).pop()
            assert local_front!=None
            
            front_verts = edge_num_to_vert_pair[local_front]
            image_verts = set(fembed.vertices().inverse()[v] for v in front_verts) # map tet verts to triangle verts
            image_edge = (set(range(3))-image_verts).pop()
            front_edges_in_internal_face.update({image_edge})
            
        if len(front_edges_in_internal_face)==1:
            pass
        else:
            return False
    return True

def is_z2_merge(tri,angle,front_edges,check_taut=False):
    ''' Determines if tri, angle, front edges determine a taut merging structure '''
    # tri should be taut
    if check_taut:
        if not is_taut(tri,angle):
            return False
    elif not is_z2_taut(tri,angle):
        return False

    # front_edges should not be where the angles are
    for i in range(len(angle)):
        if front_edges[i] in {angle[i],5-angle[i]}:
            return False

    # now check if everything matches over faces
    # every face-pairing preserved the front edge assigned locally by pis and front-edges 
    for face in tri.triangles():
        front_edges_in_internal_face=set()
        for fembed  in face.embeddings():
            face_num = fembed.face()
            tet_num = fembed.simplex().index()
            tet_fedge = front_edges[tet_num]
            ang = angle[tet_num]
            pis = [ang,5-ang]
            # need to determine if it is a pi or the front edge that endows the front, here
            face_edges = face_num_to_edges_dict[face_num]
            local_front=None
            if tet_fedge in face_edges:
                local_front = tet_fedge
            else:
                local_front = set(pis).intersection(set(face_edges)).pop()
            assert local_front!=None
            
            front_verts = edge_num_to_vert_pair[local_front]
            image_verts = set(fembed.vertices().inverse()[v] for v in front_verts) # map tet verts to triangle verts
            image_edge = (set(range(3))-image_verts).pop()
            front_edges_in_internal_face.update({image_edge})
            
        if len(front_edges_in_internal_face)==1:
            pass
        else:
            return False
    return True  

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

def unadjusted_merge_sig_from_tri_angle(tri,angle,front_edges,return_isom=True):
    ''' Gives a triangulation's isosig in its possibly non-minimal, only applying the transformation to the regina
    sig to the angle structure. Returns the isom as well. Since we need to work with isoms, we do this so we have
    less overhead when exploring (ambiguous) move routes. The truncated taut sig still works as intended with isosig_to_tri_angle'''
    
    sig, isom = tri.isoSigDetail()
    corrected_angle = apply_isom_to_angle_struct_list(angle,isom)
    corrected_front_edges = apply_isom_to_front_edges(front_edges,isom)
    return sig+'_'+''.join([str(ang) for ang in corrected_angle])+'_'+''.join([str(fedge) for fedge in corrected_front_edges]), isom

def merge_matches_across_face(face,tri,angle,merge):
    ''' Checks if the arrows determined by a tet-front-edge assignment is well defined across the face '''
    triangle = tri.triangle(face)
    if triangle.isBoundary():
        return True
        
    front_edges = []
    for fembed in triangle.embeddings():
        tet_num =fembed.simplex().index()
        internal_face = fembed.face()
        edges = face_num_to_edges_dict[internal_face]
        front_edge = merge[tet_num]
        if front_edge in edges:
            large_edge = front_edge
        else:
            # its whatever your pi edge is
            pi_edges = edge_pair_to_edge_numbers_dict[angle[tet_num]]
            large_edge = pi_edges[0] if pi_edges[0] in edges else pi_edges[1]

        inv = fembed.vertices().inverse()
        v0, v1 = edge_num_to_vert_pair[large_edge]
        missing_vert = 3 - inv[v0] - inv[v1]
        front_edges.append(missing_vert)
    return len(set(front_edges))==1 # works with unglued faces to

def merge_matches_across_all_faces(tri,angle,merge):
    ''' check all faces for the above'''
    for face_num in range(tri.countTriangles()):
        if not merge_matches_across_face(face_num,tri,angle,merge):
            return False
    return True

def face_in_tet_num_to_front_verts(face_num,tet_num, angle, merge):
    pis = angle[tet_num],5-angle[tet_num]
    pi_verts = [set(edge_num_to_vert_pair[p]) for p in pis]
    merge_verts = set(edge_num_to_vert_pair[merge[tet_num]])
    verts = set(r for r in range(4) if r!=face_num)
    if merge_verts.issubset(verts):
        return list(merge_verts)
    else:
        return list(next(p for p in pi_verts if p.issubset(verts)))

def angle_merge_are_compatible(angle,merge):
    return all(merge[i] not in (angle[i],5-angle[i]) for i in range(len(angle)))
    
def merge_matches_taut_angle(tri,angle,merge):
    ''' checks if the merge's dual structure is the given angle'''
    face_large_edges = merge_to_face_large_edges(tri,angle,merge)
    for edge_split in edge_sides_splitted(tri,angle):
        edge_check=[]
        for face_fan in edge_split:
            fan_check = []
            for face_pair in face_fan:
                face_num,edge_num = face_pair
                coor = int(face_large_edges[face_num]==edge_num)
                fan_check.append(coor)
            if len(set(fan_check))!=1:
                # the directions dont change at pis
                return False
            edge_check.append(fan_check[0])
        if edge_check[0]==edge_check[1]:
            # it's not in and out
            return False
    return True

def apply_isom_to_front_edges(merge,isom):
    ''' Given a merge on tri, apply a regina isom to it'''
    output=[None]*len(merge)
    for i in range(len(merge)):
        tet_num = isom.tetImage(i)
        front_verts = edge_num_to_vert_pair[merge[i]]
        mapped = [isom.facetPerm(i)[front_verts[0]], isom.facetPerm(i)[front_verts[1]]]
        mapped.sort()
        mapped = tuple(mapped)
        output[tet_num] = vert_pair_to_edge_num[mapped]
    return output

def get_merge_tet_types(tri,angle,merge):
    ''' Returns a list whose ith entry is the parity of the triple point dual to tet i'''
    assert tri.isOriented() # otherwise left and right are not well-defined
    tet_types = []
    for i in range(tri.countTetrahedra()):
        tet = tri.tetrahedron(i)
        assert tet.orientation()==1
        front = merge[i]
        large_edgepair = angle[i]
        edge_pair = min(merge[i],5-merge[i])
        # now use the orientation, invert the tet, and use the veering convention:
        if large_edgepair == (edge_pair + 1) % 3:
            tet_types.append('L')
        else:
            assert large_edgepair == (edge_pair + 2) % 3
            tet_types.append('R')
    return tet_types

def is_monochromatic(tri,angle,merge):
    return len(set(get_merge_tet_types(tri,angle,merge)))==1
    
''' WARNING: the following are adjustments from veering to the z2 merge setting=========================================='''
   
def lex_smallest_angle_merge_structure(tri, angle, merge, return_isom = False):
    """
    Finds the lexicographically smallest angle structure among
    symmetries of the one we have.
    """
    all_isoms = tri.findAllIsomorphisms(tri)
    all_pairs = []
    original_angle = angle[:]
    original_merge = merge[:]
    #assert is_merge(tri,original_angle,original_merge)
    for isom in all_isoms:
        angle = apply_isom_to_angle_struct_list(original_angle, isom) 
        merge = apply_isom_to_front_edges(original_merge,isom)
        all_pairs.append( (angle,merge, isom) )
    all_pairs.sort(key = lambda pair: pair[:-1])
    if return_isom:
        return all_pairs[0]
    else:
        return all_pairs[0][0:-1]

def append_front_edges_to_sig(tri,angle,merge,return_isom=False):
    isosig, isosig_isom = tri.isoSigDetail()  # isom is the mapping from the original triangulation to the isosig triangulation
    isosig_tri = regina.Triangulation3.fromIsoSig(isosig)
    angle = apply_isom_to_angle_struct_list(angle, isosig_isom)
    merge = apply_isom_to_front_edges(merge,isosig_isom)
    min_angle,min_merge,isom2 = lex_smallest_angle_merge_structure(isosig_tri,angle,merge,return_isom=True)
    sig = tri.isoSig()+'_'+''.join(str(a) for a in min_angle)+'_'+''.join(str(m) for m in min_merge)
    if not return_isom:
        return sig
    else:
        return sig,isom2 * isosig_isom

def is_transverse_taut(triangulation, taut_angle_struct, return_type = "boolean"):
    ''' [LH]: ripped and changed from veering.transverse_taut to work for Z/2'''
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

#==== some optimisations for better merge_census-from-[blank]_census type studies
def find_all_transitive_z_merge_structures(tri, angle=None):
    if angle is not None:
        return find_transitive_z2_merge_structures_fast(tri, angle)
    results = []
    for a in get_taut_structures(tri):
        for merge in find_transitive_z2_merge_structures_fast(tri, a):
            results.append((a, merge))
    return results
        
def try_propagate(start_choice,edge_count,edge_choice_to_out_pairs,in_pair_to_forced):
    # defined here lazily... TODO: move out
    assignment = [-1] * edge_count
    queue = [(0, start_choice)] # assumes first coor is zero
    
    while queue:
        e, choice = queue.pop()
        if assignment[e] != -1:
            if assignment[e] == choice:
                continue
            else:
                return None  # contradiction
        assignment[e] = choice
        
        for out_pair in edge_choice_to_out_pairs[(e, choice)]:
            face_num, internal_edge = out_pair
            for in_edge in range(3):
                if in_edge != internal_edge:
                    in_pair = (face_num, in_edge)
                    if in_pair in in_pair_to_forced:
                        forced_edge, forced_choice = in_pair_to_forced[in_pair]
                        queue.append((forced_edge, forced_choice))
    
    if -1 in assignment:
        return None  # when we've exhausted the current sc-component
    
    return assignment

def find_transitive_z2_merge_structures_fast(tri, angle):
    # for transitive structures, we don't need any DFS stuff, so we can just fell-swoop from an initial coorientation
    # not that being transitive is not guaranteed from this, but it will not guarantee return some non-transitive structures
    # for example if you set an edge coor on a unique maximal strongly component (if it exists), this will also show up. 
    edge_count = tri.countEdges()
    edge_splits = edge_sides_splitted(tri, angle)
    both_coors = [load_both_coor_choices(edge_splits[e]) for e in range(edge_count)]
    
    in_pair_to_forced = {}
    for e in range(edge_count):
        edge_split = edge_splits[e]
        for choice_idx in range(2):
            coor = both_coors[e][choice_idx]
            in_slots = [i for i in range(len(coor)) if coor[i] == 0]
            for slot in in_slots:
                for pair in edge_split[slot]:
                    in_pair_to_forced[pair] = (e, choice_idx)
    
    edge_choice_to_out_pairs = {}
    for e in range(edge_count):
        edge_split = edge_splits[e]
        for choice_idx in range(2):
            coor = both_coors[e][choice_idx]
            out_slots = [i for i in range(len(coor)) if coor[i] == 1]
            out_pairs = []
            for slot in out_slots:
                out_pairs.extend(edge_split[slot])
            edge_choice_to_out_pairs[(e, choice_idx)] = out_pairs
    
    results = []
    for start_choice in range(2):
        assignment = try_propagate(start_choice, edge_count, edge_choice_to_out_pairs, in_pair_to_forced)
        if assignment is not None:
            coor = {e: both_coors[e][assignment[e]] for e in range(edge_count)}
            if is_z2_merge_coor(coor, tri, angle):
                merge = coor_mapping_to_front_edges(tri, angle, coor)
                results.append(merge)
    
    return results

def get_transitive_merge_sigs(tri):
    return list(set(append_front_edges_to_sig(tri,angle,merge) for angle,merge in find_all_transitive_z_merge_structures(tri)))