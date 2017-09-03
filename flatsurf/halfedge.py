import numpy as np

"""
Classes and functions related to creating HalfEdgeTriangulation data structures.
"""


class Vertex(object):
    def __init__(self, halfedge=None, coord_3d=None,
                 coord_2d_orig=None, coord_2d_prime=None):
        # Only need pointer to one halfedge
        self.halfedge = halfedge
        self.coord_3d = coord_3d
        self.coord_2d_orig = coord_2d_orig
        self.coord_2d_prime = coord_2d_orig if coord_2d_prime is None else coord_2d_prime

    @property
    def coord_2d(self):
        if self.coord_2d_prime is None:
            self.coord_2d_prime = self.coord_2d_orig
        return self.coord_2d_prime

    def get_halfedges(self):
        e_start = self.halfedge
        edges = [e_start]
        e = e_start.e_op.e_prev
        while e != e_start:
            edges.append(e)
            e = e_start.e_op.e_prev
        return edges

    def get_faces(self):
        edges = self.get_halfedges()
        faces = [e.face for e in edges]
        return faces


class Face(object):
    def __init__(self, halfedge=None, boundary=False):
        # only need pointer to a single halfedge that points to Vertex
        self.halfedge = halfedge
        self.boundary = boundary

    def get_vertices(self):
        # Returns vertices for a face.
        e_start = self.halfedge
        vs = [e_start.vertex]
        e = e_start.e_next
        while e != e_start:
            vs.append(e.vertex)
            e = e.e_next
        return vs

    def get_halfedges(self):
        e_start = self.halfedge
        edges = [e_start]
        e = e_start.e_next
        while e != e_start:
            edges.append(e)
            e = e.e_next
        return edges

    def get_adjacent_faces(self):
        # Get shared faces in available
        e_start = self.halfedge
        faces = [e_start.e_op.face]
        e = e_start.e_next
        while e != e_start:
            faces.append(e.e_op.face)
            e = e.e_next
        return faces

    def get_shared_edge(self, face):
        # Given two faces, return shared edge
        if self == face:
            return ValueError('Faces are the same!')
        # Get halfedges and check face pointer of e_op.
        es = self.get_halfedges()
        for e in es:
            if e.e_op.face == face:
                return e
        return None


class HalfEdge(object):
    def __init__(self, vertex=None, e_op=None, e_next=None,
                 e_prev=None, face=None):
        self.vertex = vertex
        self.e_op = e_op
        self.e_next = e_next
        self.e_prev = e_prev
        self.face = face


class HalfEdgeTriangulation(object):
    def __init__(self, vertices, faces, boundaries=None):
        self.vertices = vertices
        self.faces = faces
        self.boundaries = [] if boundaries is None else boundaries

    def _edges_from_simplex(self, simplex):
        """Get list of edges for simplex.

        simplex: List of three integers corresponding to self.vertices indices.
        """
        return [tuple(
            sorted(
                (self.vertices[simplex[i]],
                 self.vertices[simplex[(i + 1) % 3]])))
                for i in range(3)]

    def _create_connectivity_dicts(self):
        self.vertex_to_face_dict = {v: {'faces': [], 'halfedges': []}
                                    for v in self.vertices}
        self.face_to_vertex_dict = {}

    def _create_edge_set_and_faces_from_simplices(self, simplices):
        """Create edge set and faces from simplices and add info to dicts.

        Args:
            simplices: List of list of integer indices for vertices of a simplex.
        """
        self.edge_set = set()
        self.faces = []
        for simplex in simplices:
            # Add new edges to edge_set.
            new_edges = self._edges_from_simplex(simplex)
            self.edge_set = self.edge_set.union(set(new_edges))
            self.faces.append(Face())
            self.face_to_vertex_dict.update({self.faces[-1]: {'vertices': []}})
            for vertex_index in simplex:
                # For each vertex, append face.
                self.vertex_to_face_dict[
                    self.vertices[vertex_index]]['faces'].append(self.faces[-1])
                # For each face, append vertex.
                self.face_to_vertex_dict[
                    self.faces[-1]]['vertices'].append(self.vertices[vertex_index])

    def _create_halfedges(self):
        """For each edge in self.edge_set, create two HalfEdges and add e_op pointers.
        """
        for vertex_1, vertex_2, in self.edge_set:
            halfedge_1 = HalfEdge(vertex_1)
            halfedge_2 = HalfEdge(vertex_2)
            # Every vertex only stores one incident halfedge.
            if vertex_1.halfedge is None: vertex_1.halfedge = halfedge_1
            if vertex_2.halfedge is None: vertex_2.halfedge = halfedge_2
            halfedge_1.e_op, halfedge_2.e_op = halfedge_2, halfedge_1
            self.vertex_to_face_dict[vertex_1]['halfedges'].append(halfedge_1)
            self.vertex_to_face_dict[vertex_2]['halfedges'].append(halfedge_2)

    def _get_face_candidate_halfedges(self, face):
        """Get candidate HalfEdges of Face.

        Args:
            face: Face to get candidate halfedges for.

        Returns:
            List of HalfEdges.
        """
        vertices = self.face_to_vertex_dict[face]['vertices']
        # Get list of halfedges with vertices on face.
        halfedges = [halfedge for vertex in vertices
                     for halfedge in self.vertex_to_face_dict[vertex]['halfedges']
                     if halfedge.e_op.vertex in vertices]
        if len(halfedges) != 6:
            raise ValueError('There are %d instead of 6 halfedges' % len(halfedges))
        else:
            return halfedges

    def _set_face_pointer_if_none(self, face, halfedges):
        if face.halfedge == None:
            halfedge = halfedges[0]
            face.halfedge = halfedge

    def _get_prev(self, halfedge, halfedges):
        # Want halfedge that points to the same as opposite of halfedge
        #  but is not the opposite of h.
        for previous_halfedge in halfedges:
            if (previous_halfedge.vertex == halfedge.e_op.vertex
                and previous_halfedge != halfedge.e_op):
                return previous_halfedge
            else:
                continue

    def _set_prev_of_halfedge_and_next_of_halfedge_prev(self, halfedge, halfedges):
        previous_halfedge = self._get_prev(halfedge, halfedges)
        halfedge.e_prev, previous_halfedge.e_next = previous_halfedge, halfedge
        return halfedge.e_prev

    def _get_halfedge_candidate_faces(self, halfedge):
        faces_1 = self.vertex_to_face_dict[halfedge.vertex]['faces']
        faces_2 = self.vertex_to_face_dict[halfedge.e_op.vertex]['faces']
        candidate_faces = list(set(faces_1).intersection(set(faces_2)))
        # Candidate faces should be either 2 or 1 (if on boundary).
        if len(candidate_faces) != 2 or len(candidate_faces) != 1:
            return candidate_faces
        else:
            raise ValueError

    def _set_up_halfedge_face_if_not_on_boundary(self, halfedge, face, faces):
        other_face = faces[0] if faces[0] != face else faces[1]
        # Assign face of halfedge
        halfedge.e_op.face = other_face
        if other_face.halfedge is None:
            # Assign halfedge of face
            other_face.halfedge = halfedge.e_op
        # Add other_f to queue if in face_set and remove
        if other_face.boundary:
            self.face_queue.append(other_face)
        if other_face in self.face_set:
            self.face_queue.append(other_face)
            self.face_set.remove(other_face)

    def _set_prev_and_next_pointers(self, face, halfedges):
        halfedge = face.halfedge
        halfedge_start = halfedge
        while True:
            halfedge.face = face
            if halfedge.e_op.face is None:
                # Set up op face
                faces = self._get_halfedge_candidate_faces(halfedge)
                # If opposite half edge not on boundary
                if len(faces) == 2:
                    self._set_up_halfedge_face_if_not_on_boundary(halfedge, face, faces)
                # Handle case where other_f is boundary.
                elif len(faces) == 1:
                    self.boundary_halfedges_set.add(halfedge.e_op)
            # Set up prev of halfedge and next of halfedge prev
            halfedge = self._set_prev_of_halfedge_and_next_of_halfedge_prev(halfedge, halfedges)
            # Break loop if we return to halfedge_start
            if halfedge == halfedge_start:
                break

    def _create_boundary_faces(self):
        while len(self.boundary_halfedges_set) > 0:
            # Get a boundary halfedge.
            halfedge = self.boundary_halfedges_set.pop()
            # Make a boundary face.
            face = Face(boundary=True)
            # Add boundary face to boundaries list.
            self.boundaries.append(face)
            halfedge_start = halfedge
            halfedge.face = face
            while True:
                # Current boundary edge
                current_halfedge = halfedge
                # Find next halfedge
                while halfedge.face != None:
                    if halfedge.face == face:
                        halfedge = halfedge.e_op.e_prev.e_op
                        if halfedge.face == face:
                            break
                    else:
                        halfedge = halfedge.e_prev.e_op
                        if halfedge.face == face:
                            break
                self.boundary_halfedges_set.discard(halfedge)
                current_halfedge.e_next, halfedge.e_prev = halfedge, current_halfedge
                halfedge.face = face
                if halfedge == halfedge_start:
                    break

    def _connect_halfedges_faces_and_vertices(self):
        """Connect HalfEdges, Faces, and Vertices
        """
        self.face_set = set(self.faces)
        # Initial face is simply the first popped off self.faces.
        face = self.face_set.pop()
        self.face_queue = [face]

        self.boundary_halfedges_set = set()

        while len(self.face_queue) > 0:
            face = self.face_queue.pop(0)
            halfedges = self._get_face_candidate_halfedges(face)
            # Check if halfedges of face is not defined.
            # This should only occur for the first triangle
            self._set_face_pointer_if_none(face, halfedges)
            # Set prev and next pointers
            self._set_prev_and_next_pointers(face, halfedges)

        self._create_boundary_faces()

    @classmethod
    def create_triangulation_from_coords_and_simplices(cls, coords, simplices):
        """Create triangulation from coordinates of vertices and simplices.

        coords: List of np.array of shape (3,), coordinates for vertices.
        simplices: List of triples of indices of coords signifying which vertices
            are triangles.
        """
        triangulation = cls([], [], None)
        triangulation.vertices = [Vertex(coord_3d=coord) for coord in coords]

        # Initialize connectivity dictionaries.
        triangulation._create_connectivity_dicts()
        # Create edge set and faces from simplices
        triangulation._create_edge_set_and_faces_from_simplices(simplices)

        # Create halfedges from edge set and add to connectivity dictionaries.
        triangulation._create_halfedges()

        # Connect triangulation HalfEdges, Faces, and Vertices
        triangulation._connect_halfedges_faces_and_vertices()

        # Delete unused parameters
        delattr(triangulation, 'vertex_to_face_dict')
        delattr(triangulation, 'face_to_vertex_dict')
        delattr(triangulation, 'edge_set')
        delattr(triangulation, 'face_set')
        delattr(triangulation, 'face_queue')
        delattr(triangulation, 'boundary_halfedges_set')

        return triangulation

    def get_edges_to_vertex(self, halfedge_incident):
        halfedge_start = halfedge_incident
        current_halfedge = halfedge_start.e_next
        current_halfedge = current_halfedge.e_op
        edges = [halfedge_start]
        while current_halfedge != halfedge_start:
            edges.append(current_halfedge)
            current_halfedge = current_halfedge.e_next
            current_halfedge = current_halfedge.e_op
        return edges