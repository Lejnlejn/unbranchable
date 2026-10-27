# veering_exotic.py
# Basic functions for obtaining the exotic taut structure from a veering triangulations
from veering.taut import isosig_from_tri_angle
from veering.transverse_taut import is_transverse_taut, get_top_and_bottom_vert_nums
from veering.veering_tri import is_veering
from veering.taut import unsorted_vert_pair_to_edge_pair, vert_pair_to_edge_num, isosig_to_tri_angle
from .branchable import append_front_edges_to_sig


def exotic_taut_merge(tri,angle, parity='lower'):
    # returns the isosig of the parity (=upper or lower) exotic structure
    # this also works with `wavering triangulation'
    assert tri.isOriented() # else wtf is veering bro
    colours = is_veering(tri,angle,return_type='veering_colours')
    coors = is_transverse_taut(tri, angle, return_type ='tet_vert_coorientations')
    new_angle = []
    merge=[]
    # assign a new angle for each simplex
    for tet in tri.simplices():
        vert_nums = get_top_and_bottom_vert_nums(coors, tet.index())
        top_vert_nums = tuple(vert_nums[0])
        bot_vert_nums = tuple(vert_nums[1])
        top_edge_num = vert_pair_to_edge_num[top_vert_nums]
        bot_edge_num = vert_pair_to_edge_num[bot_vert_nums]

        # sort them into equatorial pairs
        pi_edge_pair = unsorted_vert_pair_to_edge_pair[top_vert_nums]
        eq_edges = [0,1,2]
        eq_edges.remove(pi_edge_pair)

        colour_list = []
        for i in range(6):
            colour = colours[tet.edge(i).index()]
            colour_list.append(colour)
        colour_list

        top_colour = colour_list[top_edge_num]
        bot_colour = colour_list[bot_edge_num]
        eq_colours = [colour_list[x] for x in eq_edges]
        eq_mapping = dict(zip(eq_colours, eq_edges))
        
        if parity == 'upper':
            merge.append(top_edge_num)
            swap_colour = bot_colour
            
        if parity == 'lower':
            merge.append(bot_edge_num)
            swap_colour = top_colour

        new_pi = eq_mapping[swap_colour]
        new_angle.append(new_pi)

    # now make each edge the upper or lower pi
    return new_angle, merge
        
def exotic_taut_merge_sig(sig, parity='lower', return_type='sig'):
    # returns the isosig of the parity (=upper or lower) exotic structure
    # this always works with `wavering triangulation'
    tri, angle = isosig_to_tri_angle(sig) # so it is oriented
    colours = is_veering(tri,angle,return_type='veering_colours')
    coors = is_transverse_taut(tri, angle, return_type ='tet_vert_coorientations')
    new_angle = []
    merge=[]

    for tet in tri.simplices():
        vert_nums = get_top_and_bottom_vert_nums(coors, tet.index())
        top_vert_nums = tuple(vert_nums[0])
        bot_vert_nums = tuple(vert_nums[1])
        top_edge_num = vert_pair_to_edge_num[top_vert_nums]
        bot_edge_num = vert_pair_to_edge_num[bot_vert_nums]

        # sort them into equatorial pairs
        pi_edge_pair = unsorted_vert_pair_to_edge_pair[top_vert_nums]
        eq_edges = [0,1,2]
        eq_edges.remove(pi_edge_pair)

        colour_list = []
        for i in range(6):
            colour = colours[tet.edge(i).index()]
            colour_list.append(colour)
        colour_list

        top_colour = colour_list[top_edge_num]
        bot_colour = colour_list[bot_edge_num]
        eq_colours = [colour_list[x] for x in eq_edges]
        eq_mapping = dict(zip(eq_colours, eq_edges))
        
        if parity == 'upper':
            merge.append(top_edge_num)
            swap_colour = bot_colour
            
        if parity == 'lower':
            merge.append(bot_edge_num)
            swap_colour = top_colour

        new_pi = eq_mapping[swap_colour]
        new_angle.append(new_pi)

    # now make each edge the upper or lower pi
    if return_type=='sig':
        return append_front_edges_to_sig(tri,new_angle,merge)
    else:
        return new_angle, merge

def bizzaro_exotic(sig, parity='lower'):
    # troubleshooting the above, its what happens when you confuse parities in the construction. The result is usually nonsense
    # build the tri angle
    tri, angle = isosig_to_tri_angle(sig)
    colours = is_veering(tri,angle,return_type='veering_colours')
    coors = is_transverse_taut(tri, angle, return_type ='tet_vert_coorientations')
    new_angle = []
    
    for tet in tri.simplices():
        vert_nums = get_top_and_bottom_vert_nums(coors, tet.index())
        top_vert_nums = tuple(vert_nums[0])
        bot_vert_nums = tuple(vert_nums[1])
        top_edge_num = vert_pair_to_edge_num[top_vert_nums]
        bot_edge_num = vert_pair_to_edge_num[bot_vert_nums]

        # sort them into equatorial pairs
        pi_edge_pair = unsorted_vert_pair_to_edge_pair[top_vert_nums]
        eq_edges = [0,1,2]
        eq_edges.remove(pi_edge_pair)

        colour_list = []
        for i in range(6):
            colour = colours[tet.edge(i).index()]
            colour_list.append(colour)
        colour_list

        top_colour = colour_list[top_edge_num]
        bot_colour = colour_list[bot_edge_num]
        eq_colours = [colour_list[x] for x in eq_edges]
        eq_mapping = dict(zip(eq_colours, eq_edges))
        
        if parity == 'lower':
            anti_colour = {'red','blue'}-{bot_colour}
            swap_colour = anti_colour.pop()
            
        if parity == 'upper':
            anti_colour = {'red','blue'}-{top_colour}
            swap_colour = anti_colour.pop()

        new_pi = eq_mapping[swap_colour]
        new_angle.append(new_pi)

    return isosig_from_tri_angle(tri,new_angle)