# For drawing taut_layered_ball objects from custom_continents.py, 
# We thank ChatGPT 5.6 sol for help in how to use matplotlib to 
# plot the computed order in our desired format. 

import math
import os

import matplotlib
matplotlib.use('PDF')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import numpy as np

from .continent_tracks import LandscapeTrack
from .cusp_order import compute_scores_via_sort, orient_ball_edges, update_order_incremental

def disk_layout(ball, parity):
    # Place the equatorial vertices on the unit circle.

    ball.build_circles(parity=parity, check=False)
    boundary = list(ball.vert_circle)
    # if 0 in boundary:
    anchor_index = boundary.index(0)
    boundary = boundary[anchor_index:] + boundary[:anchor_index]

    anchor = boundary[0]
    other_vertices = [vertex for vertex in boundary if vertex != anchor]
    if other_vertices:
        reference = min(other_vertices)

        if boundary.index(reference) > len(boundary) // 2:
            boundary = [anchor] + boundary[1:][::-1]

    positions = {}
    for index, vertex in enumerate(boundary):
        angle = 2.0 * math.pi * index / len(boundary)
        positions[vertex] = (math.cos(angle), math.sin(angle))

    return positions, boundary

def geodesic_segment(point_a, point_b, sample_count=100):
    # Points along the geodesic segment in \bbH^2 joining p and q.
    point_a = np.asarray(point_a, float)
    point_b = np.asarray(point_b, float)
    determinant = point_a[0] * point_b[1] - point_a[1] * point_b[0]

    if abs(determinant) < 1e-12: # straight line case needs care
        parameters = np.linspace(0, 1, sample_count)
        points = np.outer(1 - parameters, point_a) + np.outer(parameters, point_b)
        return points[:, 0], points[:, 1]
    
    right_a = 0.5 * (point_a.dot(point_a) + 1.0)
    right_b = 0.5 * (point_b.dot(point_b) + 1.0)
    centre = np.array([(right_a * point_b[1] - point_a[1] * right_b) / determinant, (point_a[0] * right_b - right_a * point_b[0]) / determinant])
    radius = math.sqrt(max(centre.dot(centre) - 1.0, 0.0))
    angle_a = math.atan2(point_a[1] - centre[1], point_a[0] - centre[0])
    angle_b = math.atan2(point_b[1] - centre[1], point_b[0] - centre[0])

    if angle_b - angle_a > math.pi:
        angle_b -= 2 * math.pi

    elif angle_a - angle_b > math.pi:
        angle_b += 2 * math.pi

    angles = np.linspace(angle_a, angle_b, sample_count)

    return centre[0] + radius * np.cos(angles), centre[1] + radius * np.sin(angles)

def geodesic_midpoint_and_normal(point_a, point_b):
    # Geodesic midpoint and a normal there
    x_values, y_values = geodesic_segment(point_a, point_b, sample_count=5)
    midpoint = np.array([x_values[2], y_values[2]])
    tangent = np.array([x_values[3] - x_values[1], y_values[3] - y_values[1]])
    tangent_norm = np.linalg.norm(tangent)
    assert tangent_norm > 1e-12  # distinct edge endpoints
    tangent = tangent / tangent_norm

    return midpoint, np.array([-tangent[1], tangent[0]])

def bezier_points(point_a, control_a, control_b, point_b, sample_count=32):
    #Sample points on a cubic Bezier curve
    parameters = np.linspace(0, 1, sample_count)[:, None]
    point_a = np.asarray(point_a, float)
    point_b = np.asarray(point_b, float)
    control_a = np.asarray(control_a, float)
    control_b = np.asarray(control_b, float)

    return (1 - parameters) ** 3 * point_a + 3 * (1 - parameters) ** 2 * parameters * control_a + 3 * (1 - parameters) * parameters ** 2 * control_b + parameters ** 3 * point_b

def draw_coastal_arcs(ax, ball, parity, positions, colour='tab:green'):
    # Draw source-to-sink arcs of the (parity) track wher parity=upper, lower. 
    track = LandscapeTrack(ball, parity)
    for source_edge, sink_edge in track.coastal_track_arcs():
        source = ball.tri.edge(source_edge)
        sink = ball.tri.edge(sink_edge)
        source_midpoint = geodesic_midpoint_and_normal(positions[source.vertex(0).index()], positions[source.vertex(1).index()])[0]
        sink_midpoint = geodesic_midpoint_and_normal(positions[sink.vertex(0).index()], positions[sink.vertex(1).index()])[0]
        x_values, y_values = geodesic_segment(source_midpoint, sink_midpoint)
        ax.plot(x_values, y_values, color=colour, linewidth=1.2, alpha=0.9, zorder=3)


def draw_rank_annulus(ax, boundary, vertex_ranks, rank_cmap, band_width):
    # Draw rank colours as fixed-width slivers outside the circle
    colour_map = cm.get_cmap(rank_cmap)
    max_rank = max(vertex_ranks.values()) if vertex_ranks else 0
    step = 360.0 / len(boundary)
    half_width = min(step, 1.5) / 2.0
    for index, vertex in enumerate(boundary):
        if vertex not in vertex_ranks:
            continue
        fraction = vertex_ranks[vertex] / max_rank if max_rank else 0.0
        colour = colour_map(fraction)
        centre = index * step
        ax.add_patch(Wedge((0, 0), 1.0 + band_width, centre - half_width, centre + half_width, width=band_width, facecolor=colour, edgecolor=colour, linewidth=0.5, antialiased=False, zorder=2))


def add_rank_legend(fig, rank_cmap):
    # Add the rank cmap as a legend. We add an arow in post
    colour_axis = fig.add_axes([0.30, 0.05, 0.40, 0.03])
    colour_axis.set_frame_on(False)
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    colour_axis.imshow(gradient, aspect='auto', cmap=cm.get_cmap(rank_cmap), extent=[0, 1, 0, 1])
    colour_axis.set_xticks([])
    colour_axis.set_yticks([])

# now the function drawing the landscape
def draw_landscape_pdf(ball, path, which='up', order_colour='stable', vertex_ranks=None, band_width=0.18, show_legend=True, mark_pi=True, mark_ref=True, arcs=False, figsize=None):
    # other cmaps: 'magma' for the unstable
    if order_colour == 'stable':
        rank_cmap='viridis'
    else:
        rank_cmap='magma'

    # more monotonic options are (BuGn,RdPu) for (stable,unstable)
    
    positions, boundary = disk_layout(ball, which)
    fig, ax = plt.subplots(1, 1, figsize=figsize or (6, 6))
    ax.add_patch(plt.Circle((0, 0), 1.0, fill=False, color='0.8', linewidth=0.6))
    vertex_count = len(positions)

    if vertex_ranks is not None:
        draw_rank_annulus(ax, boundary, vertex_ranks, rank_cmap, band_width)

    else:
        vertices = list(positions)
        x_values = [positions[vertex][0] for vertex in vertices]
        y_values = [positions[vertex][1] for vertex in vertices]
        marker_size = max(1.0, min(36.0, 4000.0 / max(vertex_count, 1)))
        ax.scatter(x_values, y_values, s=marker_size, c='black', linewidths=0, zorder=2)

    if mark_ref:
        # Mark the four original vertices of the base tet
        reference_tet = ball.tri.tetrahedron(0)
        reference_vertices = [reference_tet.vertex(i).index() for i in range(4)]

        for vertex in reference_vertices:
            if vertex in positions:
                x, y = positions[vertex]
                ax.plot([x], [y], marker='o', markersize=10, markerfacecolor='white', markeredgecolor='black', markeredgewidth=1.4, fillstyle='full', zorder=7)

    if arcs:
        if which=='up':
            colour = 'tab:green'
        else:
            colour = 'tab:purple'
        draw_coastal_arcs(ax, ball, which, positions,colour=colour)
        
    limit = 1.05 + band_width
    ax.set_aspect('equal')
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.axis('off')
    
    if show_legend and vertex_ranks is not None:
        add_rank_legend(fig, rank_cmap)
        
    fig.savefig(path, format='pdf')
    plt.close(fig)
    
    return path

# Finally, the systematic layer->draw for the sake of obtaining something like an animation
def run_order_movie(ball, rounds=2, generations_per_round=6, order_colour='stable',start='down', out_dir='order_figures', which='up', max_depth=None, rank_cmap='viridis', band_width=0.18, show_legend=True, mark_pi=True, mark_ref=True, arcs=False, assume_sort=True):
    os.makedirs(out_dir, exist_ok=True)
    # Run a sequence of layering everywhere and compute the corresponding orders. Each round, layer everywhere by generations_per_round.
    # After each round, compute the order and draw it on a circle.
    parities = [start, 'up' if start == 'down' else 'down']
    records = []
    vertex_order = None
    generation_count = 0
    for round_index in range(rounds):
        parity = parities[round_index % 2]

        for generation_index in range(generations_per_round):
            # perform a sequence of layerings
            generation_count += 1
            old_vertex_count = ball.tri.countVertices()
            ball.layer_in_direction(parity, 1)
            ball.convexify(parity='up')
            ball.convexify(parity='down')
            ball.build_circles(parity=which, check=False)
            ball.compute_edge_colours()
            edge_orientations = orient_ball_edges(ball)

            if assume_sort:
                # compute score fast: take transitive and symmetric closures
                scores, ranks = compute_scores_via_sort(ball, edge_orientations)
                vertex_order = None
            else:
                # compute score slow: check every pair. Only use for sanity checks. 
                scores, ranks, vertex_order = update_order_incremental(ball, edge_orientations, vertex_order, old_vertex_count)

            figure_path = os.path.join(out_dir, 'order_gen_%02d.pdf' % generation_count)
            draw_landscape_pdf(ball, figure_path, which=which, vertex_ranks=ranks, order_colour=order_colour, band_width=band_width, show_legend=show_legend, mark_pi=mark_pi, mark_ref=mark_ref, arcs=arcs)
            vertex_count = ball.tri.countVertices()
            tetrahedron_count = ball.tri.countTetrahedra()

            print('generation %2d (%s): %d tets, %d vertices -> %s' % (generation_count, parity, tetrahedron_count, vertex_count, figure_path))

            records.append({'generation': generation_count, 'parity': parity, 'tets': tetrahedron_count, 'vertices': vertex_count, 'scores': scores, 'ranks': ranks, 'figure': figure_path})

    return records