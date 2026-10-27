from veering.taut import edge_num_to_vert_pair, vert_pair_to_edge_num
from itertools import combinations

edge_pair_to_vert_pairs = {0:((0,1),(2,3)),1:((0,2),(1,3)),2:((0,3),(1,2))}

def edge_num_to_adj_faces(edge_num):
    ''' for navigating with reginas conventions'''
    verts  = edge_num_to_vert_pair[edge_num]
    other_verts = [v for v in range(4) if v not in verts] # these are also just the adjacent faces
    return other_verts

def get_fembed(tri,pair):
    tet_num,face_num = pair
    return next(f for f in tri.simplex(tet_num).triangle(face_num).embeddings() if f.face()==face_num and f.simplex().index()==tet_num)

def face_num_to_edges(face_num):
    verts = [i for i in range(4) if i!=face_num]
    return [vert_pair_to_edge_num[pair] for pair in combinations(verts,2)]

def plainsig_to_veering(sig):
    # for quick look up 
    for tautsig in get_all_tauts_quicker(sig):
        if is_veering(tautsig):
            return tautsig
    else:
        return None

def has_isolated_vertex(tri,angle=None,vert_num=0):
    for e in tri.edges():
        if [e.vertex(j).index() for j in range(2)]==[vert_num,vert_num]:
            return False
    return True

def snappy_name(M):
    # outputs this as a "more obvious to reader" single string, rather than the whole shebang
    cusps= M.num_cusps()
    if len(M.identify())>0:
        return str(M.identify().pop())[:-5*cusps]
    else:
        return None

def cut_along_face(face,tri):
    if face.isBoundary():
        return None
    fembed = face.embedding(0)
    tet = fembed.simplex()
    tet.unjoin(fembed.face())

sig_order = key = lambda s: (len(s), s) # needed to actually compare sigs appropriately...