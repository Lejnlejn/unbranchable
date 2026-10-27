# A custom version of veering's continents.py, The underlying machinery is about FFS'19's continents.

# Here, our continents only assume that the triangulation has the following propery, generalising veering. 
# Let T be taut. Call T _combed_ if every edge has a non-empty fan of tetrahedra on each side.

# An edge in a combed triangulation is thus a taut sutured ball. So, the local picture of layering over edges
# Is somewhere between veering and an arbitrary taut triangulation. 



import regina
from veering.transverse_taut import is_transverse_taut, edge_side_face_collections, get_tet_top_vert_nums
from veering.veering_tri import is_veering

from .branchable import get_fembed, is_z2_taut, is_transverse_taut as is_z2_transverse_taut
from .continent_helpers import edge_sides_w_diags, embedge_to_tet_num_pair, tet_pair_at_edge, face_num_meeting_edge, find_paired_face, get_top_pair, build_equatorial_circle, equatorial_boundary_edges, assert_combed, circle_agrees_with_reference
from .edge_orientations import face_coorientation

class lifted_tet:
    # A tet on the universal cover; shares its skeleta indices with the its "home" tet downstairs
    def __init__(self, tet, base_tet, merge2, angle2):
        self.tet = tet
        self.home = base_tet
        self.home_index = base_tet.index()
        self.merge2 = merge2
        self.angle2 = angle2

class taut_layered_ball:
    # This is our "custom continent" class. We did not directly fork off of Veering's continent.py, since the input
    # is more general than a veering triangulation. But the framework is similar. 

    def __init__(self, tri, angle, angle2, merge2):
        self.base_tri = tri
        self.base_angle = angle
        self.base_angle2 = angle2
        self.base_merge2 = merge2

        assert is_z2_taut(tri, angle2)  # secondary angle, used for edge coorientations, is at least Z2 taut
        secondary_coors = is_z2_transverse_taut(tri, angle2, return_type='tet_vert_coorientations')
        self.base_coors2 = None if secondary_coors is False else secondary_coors
        self.base_coors = is_transverse_taut(tri, angle, return_type='tet_vert_coorientations')
        assert_combed(tri, angle)
        self.edge_sides, self.tet_sides = edge_side_face_collections(self.base_tri, self.base_angle, return_tets=True, tet_vert_coorientations=self.base_coors)
        self.all_tet_sides = edge_sides_w_diags(self.base_tri, self.base_angle, self.base_coors)
        self.pis = [[e[1][0], e[3][0]] for e in self.all_tet_sides]
        n_base_edges = self.base_tri.countEdges()
        self.fan_above = [tet_pair_at_edge(self.base_tri, i, self.base_coors) for i in range(n_base_edges)]
        self.fan_below = [tet_pair_at_edge(self.base_tri, i, self.base_coors, below=True) for i in range(n_base_edges)]
        self.tri = regina.Triangulation3()
        self.lifteds = []
        self.angle = []
        self.angle2 = [] # for our applications, we add auxillary slots for a second angle and its associated merge, and transver structures
        self.merge2 = []
        self.coors = []
        self.coors2 = []
        self.ref_triple = None
        self.triangle_ages = dict()
        self.epoch = 0 # the counter of the layering procedure
        self.covered_triangles = []
        self.num_faces_flipped = 0
        self.coastal_flippable_edges = []
        self.flippable_edges = []
        self.coastal_edges = []
        self.upper_edges = []
        self.edge_colours = None
        self.geodesic_data = None
        # vert, edge and tet circle are recordings of the landscape as a circle of vertices, edges and tet-edge pairs
        self.vert_circle = None
        self.edge_circle = None
        self.tet_circle = None
        self.generation = 0
        self.quotient_cache = {}
        self.quotient_map_cache = {}
        self.ledge_proj_cache = {}

    def initiate(self, tet_num):
        # Add a lifted copy of the tet downstairs with index tet_num to the ball
        tet = self.tri.newTetrahedron()
        self.progress_generation()
        downstairs = self.base_tri.simplex(tet_num)
        ltet = lifted_tet(tet, downstairs, merge2=self.base_merge2[tet_num], angle2=self.base_angle2[tet_num])
        self.lifteds.append(ltet)
        self.angle.append(self.base_angle[tet_num])
        self.angle2.append(self.base_angle2[tet_num])
        self.merge2.append(self.base_merge2[tet_num])
        self.coors.append(self.base_coors[tet_num])
        self.coors2.append(self.base_coors2[tet_num] if self.base_coors2 is not None else None)

        if self.tri.countTetrahedra() == 1:
            face_num = next(f for f in range(4) if self.coors[0][f] == +1)
            face = tet.triangle(face_num)
            ref_triple = tuple(face.vertex(i).index() for i in range(3))
            if face_coorientation(face, self.coors) == -1:
                ref_triple = ref_triple[::-1]
            self.ref_triple = ref_triple

        for face_num in range(4):
            self.triangle_ages[tet.index(), face_num] = self.epoch

        if self.tri.countTetrahedra() == 1:
            self.coastal_edges = self.find_coastal_edges()
            self.frontier = {'up': {'coastal_flip': self.outer_edges_to_extend_over('up'), 'flip': self.pi_edges_to_flip_over('up'), 'upper': self.find_upper_edges('up')}, 'down': {'coastal_flip': self.outer_edges_to_extend_over('down'), 'flip': self.pi_edges_to_flip_over('down'), 'upper': self.find_upper_edges('down')}}
            self.set_active_parity('up')
            
        return ltet

    def set_active_parity(self, parity):
        self.active_parity = parity
        active_frontier = self.frontier[parity]
        self.coastal_flippable_edges = active_frontier['coastal_flip']
        self.flippable_edges = active_frontier['flip']
        self.upper_edges = active_frontier['upper']

    def progress_generation(self):
        # Clear stale caches
        self.generation += 1
        self.edge_colours = None
        self.geodesic_data = None
        self.quotient_cache.clear()
        self.quotient_map_cache.clear()
        self.ledge_proj_cache.clear()


    def ledge_num_to_proj_num(self, ledge_num):
        # Ball edge index -> base edge index
        cache_key = (self.generation, ledge_num)
        cached = self.ledge_proj_cache.get(cache_key)
        if cached is not None:
            return cached
        embedge = self.tri.edge(ledge_num).embedding(0)
        ltet_num, tet_edge_num = embedge_to_tet_num_pair(embedge)
        proj_tet_num = self.lifteds[ltet_num].home_index
        result = self.base_tri.simplex(proj_tet_num).edge(tet_edge_num).index()
        self.ledge_proj_cache[cache_key] = result
        return result

    def face_to_quotient(self, face_num):
        # Ball triangle index -> base triangle index.
        face = self.tri.triangle(face_num)
        fembed = face.embedding(0)
        ltet = self.lifteds[fembed.simplex().index()]
        internal_face = fembed.vertices()[3]
        return ltet.home.triangle(internal_face).index()


    def edge_sides_in_quotient(self, edge_num):
        # The fan around a ball edge as projected (base_tet, edge) pairs.
        cache_key = (self.generation, edge_num)
        cached = self.quotient_cache.get(cache_key)

        if cached is not None:
            return cached
        
        lifteds = self.lifteds
        result = []
        for embedge in self.tri.edge(edge_num).embeddings():
            result.append((lifteds[embedge.simplex().index()].home_index, embedge.face()))

        self.quotient_cache[cache_key] = result
        return result

    def edge_sides_in_quotient_mapping(self, edge_num):
        # Like edge_sides_in_quotient but as a dict mapping each projected pair back to the upstairs (ltet, edge) pair
        cache_key = (self.generation, edge_num)
        cached = self.quotient_map_cache.get(cache_key)
        if cached is not None:
            return cached
        
        lifteds = self.lifteds
        result = {}

        for embedge in self.tri.edge(edge_num).embeddings():
            t = embedge.simplex().index()
            e = embedge.face()
            result[lifteds[t].home_index, e] = (t, e)
        self.quotient_map_cache[cache_key] = result

        return result

    def edge_num_is_coastal(self, edge_num):
        # By construction will be the case if there are no pi edges
        proj_num = self.ledge_num_to_proj_num(edge_num)
        pi_fans = self.pis[proj_num]
        fan_pieces = self.edge_sides_in_quotient(edge_num)
        seen_pis = len([x for x in pi_fans if x in fan_pieces]) > 0

        return not seen_pis

    def find_coastal_edges(self):
        return [edge.index() for edge in self.tri.edges() if edge.isBoundary() and self.edge_num_is_coastal(edge.index())]
    
    def edge_num_is_upper(self, edge_num, parity='up'):
        # I.e. on the interior of the upper landscape. Lower fo parity='down'
        proj_num = self.ledge_num_to_proj_num(edge_num)
        pi_fans = self.pis[proj_num]
        fan_pieces = self.edge_sides_in_quotient(edge_num)
        num_pis = len([x for x in pi_fans if x in fan_pieces])

        if parity == 'up':
            uppermost_pi_fan = self.fan_above[proj_num]
        else:
            uppermost_pi_fan = self.fan_below[proj_num]

        has_uppermost_pi = uppermost_pi_fan in fan_pieces

        return num_pis == 1 and (not has_uppermost_pi)

    def find_upper_edges(self, parity='up'):
        return [edge.index() for edge in self.tri.edges() if edge.isBoundary() and self.edge_num_is_upper(edge.index(), parity=parity)]
    
    def layerability_from_fan(self, proj_num, fan_pieces, parity='up'):
        # Internal computation
        # The decision of how, if possible, we can layer over the edge whose fans in the ball are the input fan_pieces
        # An edge can be layered over upwards if one side has the uppermost tet in one tet fan there
        # Similarly for layering down.
        # When both uppermost (lowermost) sides are there, then the edge can be infilled over.

        pi_fans = self.pis[proj_num]
        full_fan0, full_fan1 = self.tet_sides[proj_num]
        fan_above = self.fan_above[proj_num]
        fan_below = self.fan_below[proj_num]
        any_pis = len([x for x in pi_fans if x in fan_pieces]) > 0
        
        if parity == 'up':
            if len(full_fan0) > 0:
                check0 = full_fan0[-1] in fan_pieces
            else:
                check0 = any_pis
            if len(full_fan1) > 0:
                check1 = full_fan1[-1] in fan_pieces
            else:
                check1 = any_pis
                
            both_top_fans = check0 and check1
            only_one_top_fan = (check0 or check1) and (not both_top_fans)
            
            if fan_above not in fan_pieces:
                if only_one_top_fan and (not any_pis):
                    if check0:
                        fan = full_fan0[-1] if len(full_fan0) > 0 else None
                    else:
                        assert check1  # paired face must match
                        fan = full_fan1[-1] if len(full_fan1) > 0 else None
                    return (True, 'extend', fan)
                    
                if both_top_fans:
                    if len(full_fan0) == 0:
                        return (True, 'flip', (None, full_fan1[-1]))
                    if len(full_fan1) == 0:
                        return (True, 'flip', (full_fan0[-1], None))
                    return (True, 'flip', (full_fan0[-1], full_fan1[-1]))
                    
            return (False, None, None)
            
        if parity == 'down':
            if len(full_fan0) > 0:
                check0 = full_fan0[0] in fan_pieces
            else:
                check0 = any_pis
            if len(full_fan1) > 0:
                check1 = full_fan1[0] in fan_pieces
            else:
                check1 = any_pis
                
            both_bot_fans = check0 and check1
            only_one_bot_fan = (check0 or check1) and (not both_bot_fans)
            
            if fan_below not in fan_pieces:
                if only_one_bot_fan and (not any_pis):
                    if check0:
                        fan = full_fan0[0] if len(full_fan0) > 0 else None
                    else:
                        assert check1  # paired face must match
                        fan = full_fan1[0] if len(full_fan1) > 0 else None

                    return (True, 'extend', fan)
                    
                if both_bot_fans:
                    if len(full_fan0) == 0:
                        return (True, 'flip', (None, full_fan1[0]))
                    
                    if len(full_fan1) == 0:
                        return (True, 'flip', (full_fan0[0], None))
                    
                    return (True, 'flip', (full_fan0[0], full_fan1[0]))
                    
            return (False, None, None)

    def edge_can_be_layered_over(self, edge_num, parity='up'):
            proj_num = self.ledge_num_to_proj_num(edge_num)
            fan_pieces = self.edge_sides_in_quotient(edge_num)
            return self.layerability_from_fan(proj_num, fan_pieces, parity)
    
    def edge_num_to_tet_face_data(self, edge_num, parity='up'):
        # Given an edge we will glue Get the face-pair data from the continent 
        verdict_tuple = self.edge_can_be_layered_over(edge_num, parity=parity)
        assert verdict_tuple[0]  # edge is layerable
        verdict, mode, proj_fan = verdict_tuple
        mapping = self.edge_sides_in_quotient_mapping(edge_num)
        proj_num = self.ledge_num_to_proj_num(edge_num)

        if mode == 'extend':
            # a coastal layering
            if proj_fan is None:
                proj_fan = self.fan_below[proj_num] if parity == 'up' else self.fan_above[proj_num]
            tet_num = mapping[proj_fan][0]
            proj_tet_num, internal_edge_num = proj_fan
            face_num = face_num_meeting_edge(self.base_tri, proj_tet_num, internal_edge_num, self.base_coors, top=(parity == 'up'))

            return (tet_num, face_num)
        
        if mode == 'flip':
            # an 'infill' layering
            output = []
            for fan in proj_fan:
                if fan is not None:
                    tet_num = mapping[fan][0]
                    proj_tet_num, internal_edge_num = fan
                    face_num = face_num_meeting_edge(self.base_tri, proj_tet_num, internal_edge_num, self.base_coors, top=(parity == 'up'))
                else:
                    if parity == 'up':
                        extreme_fan = self.fan_below[proj_num]
                        candidate_faces = set(range(4)) - set(get_tet_top_vert_nums(self.base_coors, extreme_fan[0]))
                    else:
                        extreme_fan = self.fan_above[proj_num]
                        candidate_faces = set(get_tet_top_vert_nums(self.base_coors, extreme_fan[0]))
                        
                    tet_num = mapping[extreme_fan][0]
                    boundary_faces = [i for i in candidate_faces if self.tri.tetrahedron(tet_num).triangle(i).isBoundary()]
                    assert len(boundary_faces) == 1  # flat side has one boundary face
                    face_num = boundary_faces.pop()

                output.append((tet_num, face_num))

            return output

    def outer_edges_to_extend_over(self, parity='up'):
        # Compute from scratch the coastal boundary edges that we can layer over
        edge_nums = []
        for edge_num in self.find_coastal_edges():
            verdict, method, fan = self.edge_can_be_layered_over(edge_num, parity=parity)
            if verdict and method == 'extend':
                edge_nums.append(edge_num)

        return edge_nums

    def pi_edges_to_flip_over(self, parity='up'):
        # Compute from scratch the infillable edges on the landscape of the given parity
        edge_nums = []
        for edge_num in self.find_upper_edges(parity=parity):
            verdict, method, fan = self.edge_can_be_layered_over(edge_num, parity=parity)
            if verdict and method == 'flip':
                edge_nums.append(edge_num)

        return edge_nums

    def glue_new_tet_over_face(self, ltet, face_num):
        # Glue one new tetrahedron across a boundary face and return the newly introduced tet
        tet0 = ltet.tet
        assert tet0.triangle(face_num).isBoundary()  # layering starts on the boundary
        proj_tri = self.base_tri
        pair0 = (ltet.home.index(), face_num)
        gluing0 = get_fembed(proj_tri, pair0).vertices()
        paired_proj = find_paired_face(proj_tri, pair0)
        fembed1_proj = get_fembed(proj_tri, paired_proj)
        gluing1 = fembed1_proj.vertices()
        proj1 = fembed1_proj.simplex()
        glue_verts = gluing1 * gluing0.inverse()

        newtet = self.initiate(proj1.index())
        tet0.join(face_num, newtet.tet, glue_verts)
        self.progress_generation()

        return newtet

    def update_new_tet_edges(self, newtet):
        # Refresh frontier data after adding a tetrahedron. We only need to look at the edges
        # that are adjacent to the tetrahedron we just added
        coastal_set = set(self.coastal_edges)
        up_upper_set = set(self.frontier['up']['upper'])
        down_upper_set = set(self.frontier['down']['upper'])
        for i in range(6):
            edge_num = newtet.edge(i).index()
            proj_num = self.ledge_num_to_proj_num(edge_num)
            fan_pieces = self.edge_sides_in_quotient(edge_num)
            pi_fans = self.pis[proj_num]
            for par in ('up', 'down'):
                verdict = self.layerability_from_fan(proj_num, fan_pieces, par)
                coastal_frontier = self.frontier[par]['coastal_flip']
                flip_frontier = self.frontier[par]['flip']
                
                if verdict[0]:
                    mode = verdict[1]
                    if mode == 'flip':
                        if edge_num not in flip_frontier:
                            flip_frontier.append(edge_num)
                        if edge_num in coastal_frontier:
                            coastal_frontier.remove(edge_num)
                            
                    else:
                        assert mode == 'extend'  # coastal moves extend
                        if edge_num not in coastal_frontier:
                            coastal_frontier.append(edge_num)
                        if edge_num in flip_frontier:
                            flip_frontier.remove(edge_num)
                else:
                    if edge_num in coastal_frontier:
                        coastal_frontier.remove(edge_num)
                    if edge_num in flip_frontier:
                        flip_frontier.remove(edge_num)
                        
            num_pis = sum((1 for x in pi_fans if x in fan_pieces))
            if num_pis == 0:
                if edge_num not in coastal_set:
                    self.coastal_edges.append(edge_num)
                    coastal_set.add(edge_num)
                    
            elif num_pis == 1:
                for side_parity, upper_set in (('up', up_upper_set), ('down', down_upper_set)):
                    uppermost = self.fan_above[proj_num] if side_parity == 'up' else self.fan_below[proj_num]
                    if uppermost not in fan_pieces and edge_num not in upper_set:
                        self.frontier[side_parity]['upper'].append(edge_num)
                        upper_set.add(edge_num)

        self.set_active_parity(self.active_parity)

    def develop_over_face(self, ltet_num, face_num, parity='up'):
        # Extends over the boundary of the continet along ltet_num, face_num
        ltet = self.lifteds[ltet_num]
        newtet = self.glue_new_tet_over_face(ltet, face_num)
        covered_face = ltet.tet.triangle(face_num)
        self.covered_triangles.append(covered_face.index())
        self.update_new_tet_edges(newtet.tet)

        return newtet

    def develop_over_pair_of_faces(self, ltet_pair0, ltet_pair1, parity='up'):
        # Infills a tetrahedra that can be infilled over ltet_pair0 and ltet_pair1 
        ltet_num0, face0 = ltet_pair0
        ltet_num1, face1 = ltet_pair1
        ltet0 = self.lifteds[ltet_num0]
        ltet1 = self.lifteds[ltet_num1]
        assert ltet0.tet.triangle(face0).isBoundary()
        assert ltet1.tet.triangle(face1).isBoundary()

        newltet = self.glue_new_tet_over_face(ltet0, face0)
        newtet = newltet.tet
        proj_tri = self.base_tri
        pair = (ltet1.home.index(), face1)
        gluing0 = get_fembed(proj_tri, pair).vertices()
        paired_proj = find_paired_face(proj_tri, pair)
        gluing1 = get_fembed(proj_tri, paired_proj).vertices()
        glue_verts = gluing1 * gluing0.inverse()
        ltet1.tet.join(face1, newtet, glue_verts)

        self.progress_generation()
        self.num_faces_flipped += 1
        self.update_new_tet_edges(newtet)

        return newtet

    def fully_layer(self, parity='up'):
        # In no particular order, fix the set of current coastal edges and infill edges and layer over all of them. 
        self.set_active_parity(parity)
        edge_nums = list(self.coastal_flippable_edges)
        flip_nums = list(self.flippable_edges)
        self.epoch += 1
        other = 'down' if parity == 'up' else 'up'
        opposite_frontier = self.frontier[other]
        
        for edge_num in edge_nums:
            self.coastal_flippable_edges.remove(edge_num)
            if edge_num in opposite_frontier['coastal_flip']:
                opposite_frontier['coastal_flip'].remove(edge_num)
                
            if edge_num in opposite_frontier['flip']:
                opposite_frontier['flip'].remove(edge_num)
                
            if edge_num in self.coastal_edges:
                self.coastal_edges.remove(edge_num)
                
            tet, face = self.edge_num_to_tet_face_data(edge_num, parity=parity)
            self.develop_over_face(tet, face, parity=parity)
            
        for edge_num in flip_nums:
            self.flippable_edges.remove(edge_num)
            if edge_num in opposite_frontier['flip']:
                opposite_frontier['flip'].remove(edge_num)
                
            if edge_num in opposite_frontier['coastal_flip']:
                opposite_frontier['coastal_flip'].remove(edge_num)
                
            if edge_num in self.upper_edges:
                self.upper_edges.remove(edge_num)
                
            if edge_num in opposite_frontier['upper']:
                opposite_frontier['upper'].remove(edge_num)
                
            pair0, pair1 = self.edge_num_to_tet_face_data(edge_num, parity=parity)
            self.develop_over_pair_of_faces(pair0, pair1, parity=parity)

    def layer_in_direction(self, parity, num_generations):
        # Run fully layer for a fixed number of generations in the direct given by parity
        self.set_active_parity(parity)
        gens = 0

        for generation_index in range(num_generations):
            if not (self.coastal_flippable_edges or self.flippable_edges):
                break
            self.fully_layer(parity=parity)
            gens += 1

        return gens

    def convexify(self, parity='up'):
        # Convexifying a continent means to exhaust all infill tetrahedra until none are available.
        # FSS19 show that this always terminates for a veering triangulations. One can show that it always terminates for 
        # an arbitrary continent. However, the bound on how many infills it is much worse.

        self.set_active_parity(parity)
        flip_set = self.flippable_edges
        other = 'down' if parity == 'up' else 'up'
        opposite_frontier = self.frontier[other]
        while flip_set:
            edge_num = flip_set.pop()
            verdict, method, fan = self.edge_can_be_layered_over(edge_num, parity=parity)
            if not (verdict and method == 'flip'):
                continue
            
            if edge_num in opposite_frontier['flip']:
                opposite_frontier['flip'].remove(edge_num)

            if edge_num in self.upper_edges:
                self.upper_edges.remove(edge_num)

            if edge_num in opposite_frontier['upper']:
                opposite_frontier['upper'].remove(edge_num)

            pair0, pair1 = self.edge_num_to_tet_face_data(edge_num, parity=parity)
            self.develop_over_pair_of_faces(pair0, pair1, parity=parity)

    def get_upper_faces(self, parity='up'):
        # Build from scratch the upper landscape faces
        faces = []
        for face in self.tri.triangles():
            if not face.isBoundary():
                continue
            is_upper = get_top_pair(face, self.coors) is None
            if parity == 'up' and is_upper or (parity != 'up' and (not is_upper)):
                faces.append(face.index())

        return faces

    # TODO: HERE
    def build_circles(self, parity='up', check=True):
        # Build tet_circle, edge_circle and vert_circle. The equatorial circle data of the continent.
        if check:
            up_edges = equatorial_boundary_edges(self, parity='up')
            down_edges = equatorial_boundary_edges(self, parity='down')
            assert up_edges == down_edges  # upper and lower equators agree

        self.tet_circle = build_equatorial_circle(self, parity=parity)
        self.edge_circle = [self.tri.simplex(p[0]).edge(p[1]).index() for p in self.tet_circle]
        vertex_circle = []
        n = len(self.edge_circle)

        for i in range(n):
            j = (i + 1) % n
            edge0 = self.tri.edge(self.edge_circle[i])
            edge1 = self.tri.edge(self.edge_circle[j])
            verts0 = {edge0.vertex(k).index() for k in range(2)}
            verts1 = {edge1.vertex(k).index() for k in range(2)}
            matching_vertices = verts0.intersection(verts1)
            assert len(matching_vertices) == 1  # each equator step shares one vertex
            vertex_circle.append(matching_vertices.pop())

        if not circle_agrees_with_reference(vertex_circle, self.ref_triple):
            vertex_circle = vertex_circle[::-1]

        self.vert_circle = vertex_circle

    def compute_edge_colours(self):
        # Get colours on the balls from the colours on its projection
        base_colours = is_veering(self.base_tri, self.base_angle, return_type='veering_colours')
        self.edge_colours = [base_colours[self.ledge_num_to_proj_num(edge_num)] for edge_num in range(self.tri.countEdges())]
        return self.edge_colours

    def edge_colour(self, edge_num):
        # Recompute veering colours if they are stale.
        if self.edge_colours is None:
            self.compute_edge_colours()

        return self.edge_colours[edge_num]

    # decides if a vertex is to the left or right of the pair edge_verts. Note that we do not actually make a consistent choice between
    # upper and lower here, so without more care we might get slightly different orders in the total order computation 
    def left_or_right(self, edge_verts, vertex):
        vertex_start, vertex_end = edge_verts
        circle = self.vert_circle
        circle_size =len(circle)

        start_index = circle.index(vertex_start)
        end_index = circle.index(vertex_end)
        vertex_index = circle.index(vertex)

        if (vertex_index - end_index) % circle_size < (start_index - end_index) % circle_size:
            return 'left'

        return 'right'