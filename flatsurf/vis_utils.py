"""
This file should contain shortcuts for quickly plotting 2d and 3d triangulations.
"""

import plotly
import plotly.plotly as py
from plotly.graph_objs import *
import os
import numpy as np

import flatsurf.halfedge as halfedge
import flatsurf.mccartney1999 as mc1999

plotly.tools.set_credentials_file(username=os.environ['plotly_user'], api_key=os.environ['plotly_key'])

DEFAULT_AXIS = dict(showbackground=True,
                    backgroundcolor="rgb(230, 230,230)",
                    gridcolor="rgb(255, 255, 255)",
                    zerolinecolor="rgb(255, 255, 255)",
                    )

DEFAULT_LAYOUT = Layout(show_ledgend=False,
                        width=500,
                        height=500,
                        scene=Scene(xaxis=XAxis(DEFAULT_AXIS),
                                    yaxis=YAxis(DEFAULT_AXIS),
                                    zaxis=ZAxis(DEFAULT_AXIS),
                                    aspectratio=dict(x=1,y=1,z=0.5),))

def edges_from_simplex(simplex):
    """Get list of edges for simplex.

    simplex: List of three integers corresponding to self.vertices indices.
    """
    return [tuple(
        sorted(
            (simplex[i],
             simplex[(i + 1) % 3])))
            for i in range(3)]


def edge_set_from_simplices(simplices):
    edge_set = set()
    for simplex in simplices:
        new_edges = edges_from_simplex(simplex)
        edge_set = edge_set.union(new_edges)
    return edge_set


def plotly_trisurf_3d(coords, simplices, edge_set=None):
    """Make 3d plot of triangulated surface.

    Args:
        coords: List of 3d coordinates.
        simplices: List of triplets of coordinate indices.
    """
    coords = np.array(coords)
    simplices = np.array(simplices)
    points = Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                       mode='text',
                       text=map(str, range(coords.shape[0])))
    triangles = Mesh3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                       i=simplices[:, 0], j=simplices[:, 1], k=simplices[:, 2])
    edge_set = edge_set_from_simplices(simplices) if edge_set is None else edge_set
    lines = []
    for vertex1, vertex2 in edge_set:
        x, y, z = zip(coords[vertex1].tolist(), coords[vertex2].tolist())
        lines += [Scatter3d(x=x, y=y, z=z, mode='lines', line=Line(color='rgb(50,50,50)', width=1.5))]
    return Data([points, triangles] + lines)