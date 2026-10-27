"""
branch_dfs.py

A tet-by-tet DFS for branchings. This is optimised a bit to run at a reasonable rate on Shana Li's massive censuses.
For more detail on what constitutes a z2_merge structure, see branchable.py
"""

import regina
import snappy
from veering.taut import edge_num_to_vert_pair, vert_pair_to_edge_num
from .branchable import is_z2_taut, is_taut, is_transverse_taut, is_z2_merge

# precomputation
tet_edges_in_face = {f: tuple(e for e in range(6) if f not in edge_num_to_vert_pair[e])
                 for f in range(4)}

edge_num_to_adj_faces = {}
for e in range(6):
    a, b = edge_num_to_vert_pair[e]
    edge_num_to_adj_faces[e] = tuple(v for v in range(4) if v not in (a, b))

def face_large_edge_in_tet(face_num, angle, merge):
    fe = tet_edges_in_face[face_num]
    if merge in fe:
        return merge
    # if the front edge is not the front, the pi edge is
    pis = (angle,5-angle)
    for p in pis:
        if p in fe:
            return p
        
    assert False # do not get here :)

def image_edge_of_tet_edge(fembed, tet_edge):
    inv = fembed.vertices().inverse()
    a, b = edge_num_to_vert_pair[tet_edge]
    ia, ib = inv[a], inv[b]
    return 3 - ia - ib  # the missing vertex of the triangle

# we only use this once...
def tet_edge_from_image_edge(fembed, image_edge):
    # given a triangle embeddeing and an edge in the triangle, give the edge in the tetrahedron
    verts = fembed.vertices()
    tri_verts = [v for v in range(3) if v != image_edge]
    a, b = verts[tri_verts[0]], verts[tri_verts[1]]
    return vert_pair_to_edge_num[(a, b)]

def edge_num_to_pi_pair(e):
    return min(e,5-e)

# now the dfs
class BranchDead(BaseException):
    """To exit out of a dead branch"""

class TetState:
    """Stores knowledge about one tetrahedron in a partially filled out merge structure"""

    def __init__(self):
        self.possible_pis = {0, 1, 2}      # pi pairs still allowed
        self.merge = None                  # known front edge of tet
        self.face_front = [None, None, None, None]  # tet front edge per face, or None
        self.assigned = False              # have we written its 4 faces yet

    def angle(self):
        if len(self.possible_pis) == 1: # then we have a pi
            return next(iter(self.possible_pis))
        return None

    def is_determined(self):
        # do we know the total angle and the merge
        return self.angle() is not None and self.merge is not None

    def set_face_front(self, face_num, tet_edge):
        prev = self.face_front[face_num]
        if prev is not None:
            if prev != tet_edge:
                raise BranchDead() # front edge disagreement, exit
            
            return False  # nothing new
        
        self.face_front[face_num] = tet_edge
        return True # something new

    def set_pi(self, pi_pair):
        # try and declare the pi
        if pi_pair not in self.possible_pis:
            raise BranchDead() 
        
        if self.possible_pis != {pi_pair}:
            self.possible_pis = {pi_pair}

        self.check_merge_vs_pi()

    def exclude_pi(self, pi_pair):
        # remove pi from candidates
        if pi_pair in self.possible_pis:
            self.possible_pis.discard(pi_pair)
            if not self.possible_pis:
                raise BranchDead() # no pis left
            
        self.check_merge_vs_pi()

    def set_merge(self, e):
        # try and set front edge of tet
        if self.merge is not None:
            if self.merge != e:
                raise BranchDead() # disagreement
            return
        
        self.merge = e
        self.check_merge_vs_pi()

    def check_merge_vs_pi(self):
        # merge edge can never be a pi edge of the (forced) angle
        if self.merge is not None and len(self.possible_pis) == 1:
            ang = next(iter(self.possible_pis))
            if self.merge in (ang,5-ang):
                raise BranchDead() #merge = pi
            
        # if merge known, it also forbids the pi pair that merge belongs to. So remove it
        if self.merge is not None:
            bad = edge_num_to_pi_pair(self.merge)
            if bad in self.possible_pis:
                self.possible_pis.discard(bad)
                if not self.possible_pis:
                    raise BranchDead() # pis empty

def apply_edge_verdict(ts, e):
    """Given a TetState ts and an internal tet-edge e whose two bounding faces
    are both known, apply the front/front, front/rear, rear/rear verdict."""
    # If you know the front in two faces, then you know:
    # - if both faces agree, that this is the front of the tet
    # - if both disagree, then if their shared edge is not a pi. We do not know who is the pi and who is the front, though
    # - if one follows the other, then we know that the one ahead gives the front edge, the one behind gives the pis. See todo
    fA, fB = edge_num_to_adj_faces[e]
    frontA = ts.face_front[fA]
    frontB = ts.face_front[fB]

    if frontA is None or frontB is None:
        return # we cannot force anything
    
    a_is_front = (frontA == e)
    b_is_front = (frontB == e)
    if a_is_front and b_is_front:
        ts.set_merge(e)                 # front-front: e is the merge edge
    elif a_is_front or b_is_front:
        ts.set_pi(edge_num_to_pi_pair(e)) # front-rear: e is a pi edge # TODO: update so that this also adds front edge
    else:
        ts.exclude_pi(edge_num_to_pi_pair(e))  # rear-rear: e not a pi

class MergeSearch:
    def __init__(self, tri):
        self.tri = tri
        self.num_tets = tri.countTetrahedra()
        self.embed = {}        # (t,f) -> fembed in tet t with face f
        self.neighbour = {}    # (t,f) -> (t2,f2) or None (bdry face)
        for t in range(self.num_tets):
            tet = tri.tetrahedron(t)
            for f in range(4):
                tri_obj = tet.triangle(f)
                fe = None
                for emb in tri_obj.embeddings(): # find matching embed
                    if emb.simplex().index() == t and emb.face() == f:
                        fe = emb
                        break
                    
                self.embed[(t, f)] = fe
                if tri_obj.isBoundary():
                    self.neighbour[(t, f)] = None
                else:
                    adj_tet = tet.adjacentTetrahedron(f)
                    adj_gluing = tet.adjacentGluing(f)
                    self.neighbour[(t, f)] = (adj_tet.index(), adj_gluing[f])

        # For the taut prune, record pi edges
        self.tet_pi_global_edges = {}   # (t, pi_pair) -> (g1, g2)
        self.edge_is_boundary = [tri.edge(i).isBoundary() for i in range(tri.countEdges())]
        for t in range(self.num_tets):
            tet = tri.tetrahedron(t)
            for pi_pair in range(3):
                e1, e2 = (pi_pair,5-pi_pair)
                self.tet_pi_global_edges[(t, pi_pair)] = (tet.edge(e1).index(), tet.edge(e2).index())

    def candidates(self, ts):
        """All (angle, merge) consistent with this TetState's current knowledge."""
        out = []
        for angle in sorted(ts.possible_pis):
            pis = (angle,5-angle)

            for merge in range(6):
                if merge in pis:
                    continue
                if ts.merge is not None and ts.merge != merge:
                    continue
                
                # consistency with known faces
                ok = True
                for f in range(4):
                    known = ts.face_front[f]
                    if known is None:
                        continue
                    
                    if face_large_edge_in_tet(f, angle, merge) != known:
                        ok = False
                        break
                    
                if ok:
                    out.append((angle, merge))
        return out

class TautTally:
    """ For finding taut Ts. Prune any edge with more than two pis"""
    def __init__(self, search):
        self.search = search
        self.pi_count = [0] * len(search.edge_is_boundary)
        self.counted = [False] * search.num_tets   # has tet t's pis been tallied yet

    def clone(self):
        c = TautTally(self.search)
        c.search = self.search
        c.pi_count = list(self.pi_count)
        c.counted = list(self.counted)
        return c

    def record_tet_pis(self, t, ts):
        # If tet t's angle is now known and not yet counted, add its two pis and prune if any non-boundary edge exceeds two."""
        if self.counted[t]:
            return
        
        ang = ts.angle()
        if ang is None:
            return
        
        self.counted[t] = True
        g1, g2 = self.search.tet_pi_global_edges[(t, ang)]
        for g in (g1, g2):
            self.pi_count[g] += 1
            if  self.pi_count[g] > 2 and not self.search.edge_is_boundary[g]:
                raise BranchDead() # Not taut

def fire_completed_edges(ts, face_num, pending, t, tally=None):
    """After face of the tet t recieves a front edge, propagate the consequences."""

    if face_num is None:
        edge_list = range(6)
    else:
        edge_list = tet_edges_in_face[face_num]

    for e in edge_list:
        fA, fB = edge_num_to_adj_faces[e]
        if ts.face_front[fA] is not None and ts.face_front[fB] is not None: # then we know something
            apply_edge_verdict(ts, e)

    if tally is not None:
        tally.record_tet_pis(t, ts)   # a verdict may have just fixed the angle

    if ts.is_determined() and not ts.assigned:
        pending.append(t) # Add it to the tets to update to assigned

def propagate(search, states, pending, tally=None):
    """Process the pending list of determined-but-unwritten tets, writing their face
    images and propagating consequences """
    while pending:
        t = pending.pop()
        ts = states[t]
        if tally is not None:
            tally.record_tet_pis(t, ts)   # angle is known once determined

        if ts.assigned:
            continue

        if not ts.is_determined():
            continue

        angle, merge = ts.angle(), ts.merge
        ts.assigned = True


        for f in range(4):
            nb = search.neighbour[(t, f)]
            tet_edge = face_large_edge_in_tet(f, angle, merge)
            # record on this side too (keeps self-consistency / self-gluings honest)
            newly = ts.set_face_front(f, tet_edge)
            if newly:
                fire_completed_edges(ts, f, pending, t, tally)

            if nb is None: # boundary face face
                continue
            
            # image edge in the shared triangle, from this tet's embedding
            img = image_edge_of_tet_edge(search.embed[(t, f)], tet_edge)
            t2, f2 = nb
            # convert that image edge to the neighbour's tet num
            nb_tet_edge = tet_edge_from_image_edge(search.embed[(t2, f2)], img)
            ts2 = states[t2]
            newly2 = ts2.set_face_front(f2, nb_tet_edge)
            if newly2:
                fire_completed_edges(ts2, f2, pending, t2, tally)

            if tally is not None:
                tally.record_tet_pis(t2, ts2)   # its angle may have just collapsed


            if ts2.is_determined() and not ts2.assigned:
                pending.append(t2)

def check_glued_faces_consistent(search, states):
    # Check every branch we found does indeed match across faces
    for t in range(search.num_tets):
        ts = states[t]

        for f in range(4):
            nb = search.neighbour[(t, f)]

            if nb is None:
                continue
            
            if ts.face_front[f] is None:
                continue
            
            img = image_edge_of_tet_edge(search.embed[(t, f)], ts.face_front[f])
            t2, f2 = nb
            ts2 = states[t2]

            if ts2.face_front[f2] is None:
                continue
            img2 = image_edge_of_tet_edge(search.embed[(t2, f2)], ts2.face_front[f2])

            if img != img2:
                raise BranchDead()


def clone_states(states):
    # Clones states to make a new branch. 
    out = []
    for ts in states:
        c = TetState()
        c.possible_pis = set(ts.possible_pis)
        c.merge = ts.merge
        c.face_front = list(ts.face_front)
        c.assigned = ts.assigned
        out.append(c)
    return out

def solution_is_valid(search, sol_angle, sol_merge, tri, taut, transverse=False):
    # Double check a completed branch
    # this is where transversality is enforced, too
    nt = search.num_tets
    st = [TetState() for i in range(nt)]

    for t in range(nt):
        st[t].possible_pis = {sol_angle[t]}
        st[t].merge = sol_merge[t]
    try:
        propagate(search, st, list(range(nt)))
        check_glued_faces_consistent(search, st)
    except BranchDead:
        return False
    
    if taut:
        if not is_taut(tri, sol_angle):
            return False
    else:
        if not is_z2_taut(tri, sol_angle):
            return False

    if transverse and not is_transverse_taut(search.tri, sol_angle):
        # The entirety of the transverse check is here right now TODO: integrate transverse taut for speed
        return False
    return True


def find_z2_merge_structures_dfs(tri, angle=None, taut=False, transverse=False):
    """
    DFS for branchings, in the form of Z/2-taut merge structures
    """
    search = MergeSearch(tri)
    nt = search.num_tets
    solutions = []

    # initialise
    init = [TetState() for i in range(nt)]
    init_tally = TautTally(search) if taut else None
    try:
        if angle is not None: #install them across the tetrahedra
            for t in range(nt):
                init[t].set_pi(angle[t])

            if init_tally is not None:
                for t in range(nt):
                    init_tally.record_tet_pis(t, init[t])

            pending = []
            propagate(search, init, pending, init_tally)
            
    except BranchDead:
        # the supplied angle (or the triangulation itself) is inconsistent
        return []

    # finally we can actually do the dfs
    def recurse(states, tally):
        # get next undetermined tet
        t = None

        for i in range(nt):
            if not states[i].is_determined():
                t = i
                break

        if t is None: # then we are done
            sol_angle = [states[i].angle() for i in range(nt)]
            sol_merge = [states[i].merge for i in range(nt)]
            solutions.append((sol_angle, sol_merge))
            return

        # otherwise, we try and add structures
        for (a, m) in search.candidates(states[t]):
            # branch is a bad name given our context...
            branch = clone_states(states)
            branch_tally = tally.clone() if tally is not None else None
            ts = branch[t]

            try:
                ts.set_pi(a)
                ts.set_merge(m)
                if branch_tally is not None:
                    branch_tally.record_tet_pis(t, ts)

            except BranchDead:
                continue
            
            pending = [t]
            try:
                propagate(search, branch, pending, branch_tally)
            except BranchDead:
                continue
            
            recurse(branch, branch_tally)

    recurse(init, init_tally)

    # sanity-check, format
    out = []
    for sol_angle, sol_merge in solutions:
        if not solution_is_valid(search, sol_angle, sol_merge, tri, taut, transverse):
            continue
        if angle is not None:
            out.append(sol_merge)
        else:
            out.append((sol_angle, sol_merge))
    return out

def has_z2_merge_structure_dfs(tri, angle=None, taut=False, transverse=False):
    search = MergeSearch(tri)
    nt = search.num_tets
    init = [TetState() for _ in range(nt)]
    init_tally = TautTally(search) if taut else None
    try:
        if angle is not None:
            for t in range(nt):
                init[t].set_pi(angle[t])
            if init_tally is not None:
                for t in range(nt):
                    init_tally.record_tet_pis(t, init[t])
            pending = []
            propagate(search, init, pending, init_tally)
    except BranchDead:
        return False
    
    # finally we can actually do the dfs
    def recurse(states, tally):
        t = next((i for i in range(nt) if not states[i].is_determined()), None)
        if t is None:
            return solution_is_valid(search, [s.angle() for s in states],
                                    [s.merge for s in states], tri, taut, transverse)
        for (a, m) in search.candidates(states[t]):
            branch = clone_states(states)
            branch_tally = tally.clone() if tally is not None else None
            ts = branch[t]
            try:
                ts.set_pi(a)
                ts.set_merge(m)
                if branch_tally is not None:
                    branch_tally.record_tet_pis(t, ts)
                propagate(search, branch, [t], branch_tally)
            except BranchDead:
                continue
            if recurse(branch, branch_tally):
                return True
        return False

    return recurse(init, init_tally)

### Some basic searching tools. 
from veering.taut import get_taut_structures

def look_for_no_branching(tri,tries=50,taut=False,transverse=False):
    M=snappy.Manifold(tri)
    sigs = set()

    for i in range(tries):
        sigs.add(M.triangulation_isosig().split('_')[0]) # always add the default sig
        M.randomize()

    for sig in sorted(sigs,key=len):
        tri=regina.Triangulation3(sig)
        tri.orient()

        if not has_z2_merge_structure_dfs(tri,taut=taut,transverse=transverse):
            if len(get_taut_structures(tri))>0:
                # then it is essential
                return sig
        
    return False

def look_for_branching(tri,tries=50,taut=False,transverse=False):
    M=snappy.Manifold(tri)
    sigs = set()

    for i in range(tries):
        sigs.add(M.triangulation_isosig().split('_')[0]) # always add the default sig
        M.randomize()

    for sig in sorted(sigs,key=len):
        tri=regina.Triangulation3(sig)
        tri.orient()
        if has_z2_merge_structure_dfs(tri,taut=taut,transverse=transverse):
            return sig
        
    return False
