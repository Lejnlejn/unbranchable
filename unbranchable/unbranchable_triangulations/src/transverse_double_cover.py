#transverse_double_cover.py
''' Build the transverse taut double cover a taut triangulation.
I couldn't find this anywhere in Veering. Note that one can also use this to compute the edge-orientable double cover
of a veering triangulation by taking the TO cover for the exotic taut structure and inverting back to the veering. But the edge-orientable
double cover is already implemented in veering.'''

import regina
from veering.taut import is_taut
from veering.transverse_taut import is_transverse_taut 
from veering.veering_tri import edge_num_to_vert_pair

class LiftedTransverseTautTet:
    # just to store a copy of the tet in the double cover
    def __init__(self, tet, tet_downstairs, pi_pair, coors):
        self.tet = tet
        self.tet_downstairs = tet_downstairs
        self.pi_pair = pi_pair
        # given the pi_pair, there are two possible choices of coors. One lift takes each one. 
        self.coors = coors
        
def angle_to_coor_choices(pi):
    hemisphere0 = edge_num_to_vert_pair[pi]
    hemisphere1 = edge_num_to_vert_pair[5-pi]
    assignment0 = [-1 for i in range(4)]
    assignment1 = [-1 for i in range(4)]
    for v in hemisphere0:
        assignment0[v]=1
    for v in hemisphere1:
        assignment1[v]=1
    return assignment0,assignment1

def load_double_cover_tets(tri,angle):
    # install the two lifts of each tet downstairs
    double_cover = regina.Triangulation3()
    double_angle = []
    lifted_tets = []
    double_coors = []
    for i, tet in enumerate(tri.tetrahedra()):
        pi_pair = angle[i]
        coor_options = angle_to_coor_choices(pi_pair)
        for j in range(2):
            coors = coor_options[j]
            double_coors.append(coors)
            double_angle.append(pi_pair)
            new_tet = double_cover.newTetrahedron()
            lifted_tets.append(LiftedTransverseTautTet(new_tet, tet, pi_pair,coors))
            
    return double_cover, lifted_tets

def glue_double_cover_tets(tri, lifted_tets):
    # use the gluings downstairs to glue upstairs. We use coors to decide when two are glued
    for i, tet in enumerate(tri.tetrahedra()):
        for face in range(4):
            # get the downstairs neighbour
            nbr_tet = tet.adjacentTetrahedron(face)
            if nbr_tet is None:
                continue
            gluing_perm = tet.adjacentGluing(face)
            nbr_face = gluing_perm[face]  # the face of nbr_tet glued to our face
            nbr_index = nbr_tet.index()

            # load the two lifted copies of tet i
            for j in range(2):
                lifted = lifted_tets[2 * i + j]
                my_coor = lifted.coors[face]

                # get the lifted copy that has opposite coor
                for k in range(2):
                    nbr_lifted = lifted_tets[2 * nbr_index + k]
                    nbr_coor = nbr_lifted.coors[nbr_face]
                    if nbr_coor == -my_coor:
                        # glue if not already glued
                        if lifted.tet.adjacentTetrahedron(face) is None:
                            lifted.tet.join(face, nbr_lifted.tet, gluing_perm)
                        break

def build_transverse_double_cover(tri, angle,merge=None,weights=None):
    # face weights correspond to track weights downstairs
    # we will want to lift these because it means we half the dimension
    # of the linear program to find a point in the polytope of flow tracks
    if is_transverse_taut(tri,angle):
        # why are you calling me?
        return tri,angle

    double_tri, lifted_tets = load_double_cover_tets(tri, angle)
    glue_double_cover_tets(tri, lifted_tets)
    double_angle = [lt.pi_pair for lt in lifted_tets]
    assert is_taut(double_tri,double_angle)
    assert is_transverse_taut(double_tri,double_angle)


    output = [double_tri,double_angle]
    double_merge = None
    if merge!=None:
        double_merge = []
        for m in merge:
            double_merge.append(m)
            double_merge.append(m)
        output.append(double_merge)

    if weights!=None:
        weight_dict={}
        for f in tri.triangles():
            for fembed in f.embeddings():
                pair = fembed.simplex().index(),fembed.face()
                weight_dict[pair] = weights[f.index()]
            
        double_weights=[None for i in range(double_tri.countTriangles())] # these require some care
        for f in double_tri.triangles():
            fembed = f.embedding(0)
            lifted_tet_num,lifted_face = fembed.simplex().index(),fembed.face()
            downstairs_tet_num = lifted_tets[lifted_tet_num].tet_downstairs.index()
            pair = (downstairs_tet_num,lifted_face)
            double_weights[f.index()] = weight_dict[pair]
        output.append(double_weights)
    
    return output
