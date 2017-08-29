import numpy as np
from flatsurf.halfedge import HalfEdgeTriangulation, constructHalfEdgeFromSimplices


"""
Following
McCartney, J., Hinds, B. K. & Seow, B. L. Flattening of triangulationed surfaces incorporating darts and gussets.
CAD Comput. Aided Des. 31, 249-260 (1999).
"""

class Flattening(object):
    def __init__(self, triangulation_3d, Et=1e-3, delta=1e-3):
        self.triangulation_3d = triangulation_3d
        self.faces = self.triangulation_3d.faces
        self.vertices = self.triangulation_3d.vertices
        self.triangulations_2d = []
        self.available = set(self.faces)
        self.active = []
        self.flattened = set()
        self.Et = Et
        self.delta = delta

    def get_seed_face(self):
        faces = self.faces
        face_coords = [f.halfedge.vertex.coord_3d for f in faces]
        f_sort = zip(faces, face_coords)
        f_sort = f_sort.sort(key=lambda x: x[1])
        seed_face = f_sort[len(f_sort) / 2][0]
        return seed_face

    def flatten_face(self, face):
        # Get vertices for given face
        v = face.get_vertices()

        # Check that the face only has three vertices
        # We do not want to lay down boundaries
        # Need to handle boundary case still
        if len(v) != 3:
            raise ValueError("Face has more than three vertices.")
        coord_2d = [ver.coord_2d for ver in v]
        free = coord_2d.count(None)

        sides = [np.linalg.norm(v[0].coord_3d - v[1].coord_3d),
                 np.linalg.norm(v[1].coord_3d - v[2].coord_3d),
                 np.linalg.norm(v[2].coord_3d - v[0].coord_3d)]

        # Case 1: First triangle of a 2D triangulation
        if free == 3:
            # Place first coordinate at zero
            # May need to change this for laying down new piece...
            # COME BACK TO THIS
            v[0].coord_2d = np.array([0., 0.])
            v[1].coord_2d = np.array([sides[1], 0.])
            # intersection of circles at v[0] with r = side_3 and v[1] with r = side_2
            third_coord = get_circles_intersection(v[0].coord_2d, sides[2], v[1].coord_2d, sides[1])
            # Should return two solutions, taking first by convention
            # Also should error handle if something goes wrong
            v[2].coord_2d = third_coord[0]

        # Error handling -- this shouldn't happen.
        elif free == 2:
            raise ValueError("Something's wrong -- only 1 vertex flattened.")

        # Case 2: Some vertices already flattened.
        elif free == 1:
            i = coord_2d.index(None)
            third_coord = get_circles_intersection(v[(i + 1) % 3].coord_2d, sides[(i + 1) % 3],
                                                   v[i + 2 % 3].coord_2d, sides[(i + 2) % 3])
            # Should return two solutions, taking first by convention
            # Also should error handle if something goes wrong
            v[i].coord_2d_orig = third_coord[0]

        # Case 3: All vertices flattened.
        elif free == 0:
            # Laying down a constrained vertex

            # Need to pick a vertex (there will be two) that is exposed to current boundary of 2D triangulation
            adj_f = face.get_adjacent_faces()
            not_flattened = filter(lambda x: x in self.flattened, adj_f)

            # I'm not sure if this case should ever happen, but we're going to handle it anyway.
            if len(not_flattened) == 0:
                # Pick random point to average.
                i = np.random.randint(0, 3)
                third_coord = get_circles_intersection(v[(i + 1) % 3].coord_2d, sides[(i + 1) % 3],
                                                       v[i + 2 % 3].coord_2d, sides[(i + 2) % 3])
                v[i].coord_2d_prime = ( v[i].coord_2d + third_coord ) / 2.

            elif len(not_flattened) == 1:
                # Pick one of points on edge for constrained flattening.
                e = face.get_shared_edge(not_flattened.pop())
                i = v.index(e.vertex)

                third_coord = get_circles_intersection(v[(i + 1) % 3].coord_2d, sides[(i + 1) % 3],
                                                       v[i + 2 % 3].coord_2d, sides[(i + 2) % 3])
                v[i].coord_2d_prime = (v[i].coord_2d + third_coord) / 2.

            else:
                raise ValueError("This should've been an unconstrained flattening.")

            # Need to mark face as flattened before relaxing.
            self.flattened.add(face)

            self.relax()
            return None

        # Mark face as flattened and remove from available.
        self.flattened.add(face)
        return None

    def relax(self):
        # Relax all nodes in flattened triangulation.
        vertices = set([face.get_vertices() for face in self.flattened])
        for v in vertices:
            adjust_vertex(self.delta, self.Et, v)

    def flatten(self):
        seed = self.get_seed_face()
        self.active.append(seed)
        while len(self.active) != 0:
            face = self.active.pop()
            # Skip if face is a boundary.
            if face.boundary == True:
                continue
            self.flatten(face)
            shared = face.get_adjacent_faces()
            # Should this be a set to ensure no duplicates?
            add_to_active = [s for s in shared if s in self.available]
            # Add faces to active
            self.active.append(add_to_active)
            # Remove from available
            self.available = self.available - set(add_to_active)

### Energy Functions
def f_energy(vertex1, vertex2):
    # Calculate the energy of a given edge of the 2D triangulation
    v1p = vertex1.coord_2d
    v2p = vertex2.coord_2d
    v1 = vertex1.coord_2d_orig
    v2 = vertex2.coord_2d_orig
    return 0.5 * (np.dot(v1p, v2p) - np.dot(v1, v2)) ** 2 / np.dot(v1, v2)

def total_energy(vertex):
    # Calculate the total energy of a given node.
    energy = 0.
    edges = vertex.get_halfedges()
    v1 = vertex
    for e in edges:
        v2 = e.e_op.vertex
        energy += f_energy(v1, v2)
    return energy

def adjust_vertex(delta, Et, vertex):
    # Given a vertex, determine if adjusting the position of that vertex lowers the energy.
    adjust = [np.array([1, 0]), np.array([-1, 0]),
              np.array([0, 1]), np.array([0, -1])]
    coord_start = vertex.coord_2d
    energy_start = total_energy(vertex)
    dEs = []
    for a in adjust:
        vertex.coord_2d_prime = coord_start + delta * a
        dEs.append(energy_start - total_energy(vertex))
    change = max(dEs)
    if change > Et:
        vertex.coord_2d_prime = coord_start + delta * adjust[dEs.index(change)]
    else:
        vertex.coord_2d_prime = coord_start

### Utilities

# May want to see what the fastest implementation for this is in python and put in a separate utility file.
def get_circles_intersection(center1, radius1, center2, radius2):
    # https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
    d = np.linalg.norm(center1 - center2)
    if d > np.abs(radius1 + radius2) or d < np.abs(radius1 - radius2):
        # No solutions
        return None
    elif d < 1e-6 and np.abs(radius1 - radius2) < 1e-6:
        # Infinite solutions
        return "Circle"
    else:
        # Two solutions
        a = (radius1 ** 2 - radius2 ** 2 + d ** 2) / (2. * d)
        p2 = center1 + a * (center2 - center1) / d
        x3 = p2[0] + h(center2[1] - center1[1]) / d
        y3 = p2[1] - h(center2[0] - center1[0]) / d
        x3p = p2[0] - h(center2[1] - center1[1]) / d
        y3p = p2[1] + h(center2[0] - center1[0]) / d
        return (x3, y3), (x3p, y3p)