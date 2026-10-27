# branchable.py

"""
For working with branchings on triangulations. This breaks down to "Z/2Z taut merging structures",
which are `the' natural dual objects to branchings.

A branching on a tetrahedron is described by a taut tetrahedron with a front edge. I am no good at ASCII art, so here
is a description by words. The front edge is the edge dual to the one-sheeted model of the triple point. The two pi-angles
are given to the two local sectors around a triple point that have both in an inward and outward cooriented boundary.

A branching structure on a triangulation is encoded by two tuples: 
    - a list of pi pair for each tetrahedron. The resulting endowment of angles is a Z/2Z taut structure
    - a list of which edge in the tetrahedron is the front. 

The merging structure is also often considered as a coorientation on edges. For a Z/2Z taut structure, the pis
cut the tetrahedra about the edge into fans. A merging structure is a choice of coorientation (in or out) for each fan such
that the coorientation is opposite for to fans that meet a common pi. So, for an edge with angle 2pi, there is one in fan one out fan.
For an edge with angle 0, all edges are cooriented in or out, and so on. 
"""

import regina
from veering.taut import apply_isom_to_angle_struct_list, edge_num_to_vert_pair, vert_pair_to_edge_num
from .additional_tools import edge_pair_to_vert_pairs, get_fembed
from copy import deepcopy

# Local adaptations of Veering helpers used throughout this module.
from .veering_helpers import edge_pair_to_edge_numbers_dict, find_z2_taut_structures, fix_orientations, get_taut_structures, is_taut, is_transverse_taut, is_z2_taut, isosig_to_tri_angle

# Some precompute:
def face_num_to_edges(face_num):
    face_verts = list(set(range(4))-{face_num})
    vert_pairs = [(face_verts[0],face_verts[1]),(face_verts[0],face_verts[2]),(face_verts[1],face_verts[2])]
    return list(vert_pair_to_edge_num[pair] for pair in vert_pairs)

face_num_to_edges_dict = {i: face_num_to_edges(i) for i in range(4)}

### Tools for working with a Z/Z2 structure in the form or tri, angle, merge, a regina.Triangulation3, a Z2-taut structure, and a merging
### structure

def load_tri_angle_merge_from_sig(isosig, return_isom = False, apply_isom=None):
    """
    Given a taut isosig, returns an oriented regina triangulation and
    the list of angles for the taut angle structure, for the new
    labelling. 
    
    If return_isom, returns the isomorphism from the
    regina isosig triang to the oriented triangulation.

    If apply_isom = regina.isom3, applies isom to tri,angle,merge

    If the sig is only a taut isoSig, it will try to find a transitive merging structure, and otherwise will return an error. 

    The inverse operation is: append_front_edges_to_sig. 
    """
    data = isosig.split("_")
    assert len(data)==3
    isosig, angle, merge = data
    tri = regina.Triangulation3.fromIsoSig(isosig)
    angle = [int(i) for i in list(angle)]
    merge = [int(m) for m in list(merge)]
    isom = fix_orientations(tri, angle, merge, return_isom = return_isom)  # this does not alter what the angles_string should be
    assert tri.isOriented()
    #assert is_taut(tri, angle) # disable for generalised taut structures
    if apply_isom!=None:
        og_isom = apply_isom.inverse()
        tri = og_isom(tri)
        angle = apply_isom_to_angle_struct_list(angle, isom) 
        merge = apply_isom_to_front_edges(merge,isom)
  
    if return_isom:
        return tri, angle, merge, isom
    else:
        return tri, angle, merge

def load_both_coor_choices(edge_split):
    # given an `edge splitting' loads the two alternating assignments of 0,1s as a list'''
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

def face_num_to_large_edge(face,tri,angle,merge):
    # 'large edge' = front edg ein the face in the traint rack
    fembed = tri.triangle(face).embedding(0)
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
    # Assuming that tri,angle,merge is a valid merge structure, convert merge to the coors format
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
    # converts adge coors to an assignment of each (face,edge) pair as in or out of the edge
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
    # Determines if a mapping (face,edge): in or out determines a merge structure
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
                if len(all_chunks)>0:
                    updated_chunk = updated_chunk = current_chunk+all_chunks[0] # careful here: used to be all_chunks[0]+current_chunk
                    all_chunks = all_chunks[1:]
                    all_chunks.append(updated_chunk)
            
        if len(splitters)==0:
            all_chunks=[edge_side]
        edge_chunks.append(all_chunks)
    return edge_chunks

def face_pair_to_tet_in_between(tri,pair):
    # Given a pair of common faces about an edge, convert it to the tet pair (tet,edge) between them
    #this assumes that assume none of the edges have degree one, so that
    #the tet is uniquely identified as the tet with embeddings in both faces
    # for example, front edges->face large edges-> coors
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

def is_z2_merge(tri,angle,front_edges,check_taut=False):
    ''' Determines if tri, angle, front edges determine a (Z2)-taut merging structure'''
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
    # appending the merge to the sig given by "isosig_from_tri_angle" from veering for the Z/2Z-taut structure
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

### Depth-first searching for Z2-merge structures
# Compare to branch_dfs. Here, we use propagation of edge coors. There, we use propagation of front edges.
# Of course, these are not *really* different ways of searching. But this gives a good sanity check.
# This version is currently only implemented for triangulations without boundary edges.
# The alternate branch_dfs works for the with(out) boundary edge setting 

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
        return True, alive_coors, edge

def step_new_coors_out_one_from_a_single_coor(tri,new_edges,current_coor,edge_splits):
    ''' propagate a coor assignment to neighbouring faces'''
    next_edges = []
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
    next_edges = [new_edge_num]
    while len(next_edges)>0:
        verdict, current_coors, next_edges = step_new_coors_out_one_from_a_single_coor(
            tri, next_edges, current_coors, edge_splits
        )
        if not verdict:
            return False
    return current_coors

def find_merge_structures(tri,taut_over='Z2', return_on_hit=False):
    ''' Enumerate all structures of the given taut type'''
    return_output=[]

    if taut_over=='Z':
        angles = get_taut_structures(tri)
    else:
        assert taut_over=='Z2'
        angles = find_z2_taut_structures(tri)
    
    for angle in angles:
        angle_coor_pairs=[]
        edge_count=tri.countEdges()
        edge_splits = edge_sides_splitted(tri,angle)
        edge_num=0
        edge_split = edge_splits[edge_num]
        edge_coors = load_both_coor_choices(edge_split)
        alive_coors = [{edge_num:x} for x in edge_coors]
        unfinished_coors = [(x,edge_num) for x in alive_coors]
        
        while len(unfinished_coors)>0:
            coor,edge_num = unfinished_coors.pop()
            output = exhaust_new_coor(tri,coor,edge_num,edge_splits)
            if output!=False:
                coor = output # the one with edges proliferated, sure
                
                if len(coor.keys())==edge_count:
                    if is_z2_merge_coor(coor,tri,angle):
                        if return_on_hit:
                            return True
                        else:
                            angle_coor_pairs.append((angle,coor))

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