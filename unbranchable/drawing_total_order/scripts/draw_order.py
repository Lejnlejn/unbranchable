from orderable.custom_continents import taut_layered_ball
from orderable.branchable import isosig_to_tri_angle
from orderable.veering_exotic import exotic_taut_merge
from orderable.draw_landscape import run_order_movie

# sig = 'cPcbbbiht_12' # figure-eight knot complement's veering tri
# sig = 'eLMkbcdddhhhml_1221' # another edge-orientable eg. Remember that we need edge-orientable to see an order
sig = 'kLLALLQkbeddhjjijiihhnkhqqhabg_2111222110' # non-layered eg

# # Choose the lower exotic taut structure and draw the upper lamination using the stable order.
# tri, angle = isosig_to_tri_angle(sig)
# angle2, merge2 = exotic_taut_merge(tri, angle, parity='lower')
# ball = taut_layered_ball(tri, angle, angle2=angle2, merge2=merge2)
# ball.initiate(0)
# run_order_movie(ball, rounds=6, generations_per_round=2, start='down', which='up', order_colour='stable', out_dir=sig+'_order_movie_stable_order_upper_lamination', assume_sort=True, arcs=True)

# Choose the upper taut structure. The merging structure is the stable branched surface,
# so its edge orientation is its transverse order: the unstable order. Draw the lower lamination this way.
tri, angle = isosig_to_tri_angle(sig)
angle2, merge2 = exotic_taut_merge(tri, angle, parity='upper')
ball = taut_layered_ball(tri, angle, angle2=angle2, merge2=merge2)
ball.initiate(0)
run_order_movie(ball, rounds=8, generations_per_round=2, start='down', which='down', order_colour='unstable', out_dir=sig+'_order_movie_unstable_order_lower_lamination', assume_sort=True, arcs=True)