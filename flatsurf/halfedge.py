import numpy as np

"""
Contains classes and functions related to creating HalfEdgeTriangulation datastructures.
"""

class Vertex(object):
    def __init__(self, halfedge = None, coord_3d = None, 
                 coord_2d_orig = None, coord_2d_prime = None):
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
    def __init__(self, halfedge = None, boundary = False):
        # only need pointer to a single halfedge that points to Vertex
        self.halfedge = halfedge
        self.boundary = boundary

    def get_face_vertices(self):
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

    def get_shared_edge(self,face):
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
    def __init__(self, vertex = None, e_op = None, e_next = None, 
                 e_prev = None, face = None):
        self.vertex = vertex
        self.e_op = e_op
        self.e_next = e_next
        self.e_prev = e_prev
        self.face = face
        
class HalfEdgeTriangulation(object):
    def __init__(self, vertices, faces, boundaries = None):
        self.vertices = vertices
        self.faces = faces
        self.boundaries = [] if boundaries == None else boundaries
    def get_edges_to_vertex(self, halfedge_incident):
        start = halfedge_incident
        cur = start.e_next
        cur = cur.e_op
        edges = [start]
        while cur != start:
            edges.append(cur)
            cur = cur.e_next
            cur = cur.e_op
        return edges

    
def constructHalfEdgeFromSimplices(coords,simplices):
    """
    Given a list of 3d coordinates and triangulation simplices with indices 
    corresponding to the coordinate list, create HalfEdgeTriangulation.
    
    Note, does not handle periodic boundary conditions where multiple instances
    of the same vertex are in the same triangle.
    """
    
    # Create vertices and faces lists and edges set.
    vertices = []
    faces = []
    edges = set()
    for c in coords:
        v = Vertex(coord_3d = c)
        vertices.append(v)

    # Initialize vertex -> faces and edges dictionary.
    v2f = {v: {'faces': [] , 'halfedges' : []} for v in vertices}
    # Initialize faces -> vertices dictionary.
    f2v = {}
    for i,s in enumerate(simplices):    
        # Create face for every simplex.
        f = Face()
        faces.append(f)
        # Add to edge set.
        es = [tuple(sorted((vertices[s[k]], vertices[s[(k+1)%3]]),
                           key = lambda x : list(x.coord_3d))) for k in range(3)]
        edges = edges.union(set(es))
        # Create entry in face -> vertices dictionary for every simplex.
        f2v.update({f: {'vertices': []}})
        for j in s:
            # For each vertex, append face.
            v2f[vertices[j]]['faces'].append(f)
            # For face, append each vertex.
            f2v[f]['vertices'].append(vertices[j])
    
    triangulation = HalfEdgeTriangulation(vertices,faces)
    
    # Create HalfEdges.
    # Reference opposite edges.
    for v1,v2 in edges:
        h1 = HalfEdge(v1)
        h2 = HalfEdge(v2)
        # Ever vertex only needs to store one incident halfedge.
        if v1.halfedge == None: v1.halfedge = h1
        if v2.halfedge == None: v2.halfedge = h2
        h1.e_op, h2.e_op = h2, h1
        # Append halfedge to 
        v2f[v1]['halfedges'].append(h1)
        v2f[v2]['halfedges'].append(h2)
    
    def get_candidate_halfedges(f):
        # Given a face returns the candidate halfedges for that face
        # Should be 6 halfedges
        vs = f2v[f]['vertices']
        halfedges = [i for v in vs for i in v2f[v]['halfedges'] if i.e_op.vertex in vs]
        if len(halfedges) != 6:
            raise ValueError
        else:
            return halfedges
    
    def get_candidate_faces(h):
        # Given a halfedge returns faces adjacent to halfedge and its opposite
        f1 = v2f[h.vertex]['faces']
        f2 = v2f[h.e_op.vertex]['faces']
        candidate_faces = list(set(f1).intersection(set(f2)))
        # Should be either 2 or 1 (if on boundary)
        if len(candidate_faces) != 2 or len(candidate_faces) != 1:
            return candidate_faces
        else:
            raise ValueError
    
    def get_prev(h,hs):
        # Want halfedge that points to the same as opposite of h but is not the opposite of h.
        for p in hs:
            if p.vertex == h.e_op.vertex and p != h.e_op:
                return p
            else:
                continue
    
    # Iterate through Faces to set up next and previous edges
    face_set = set(faces)
    f = face_set.pop()
    face_queue = [f]
    
    # Initialize set of boundary halfedges
    boundary_halfedges_set = set()
    # First triangle is special because it choose the inner edge orientation convention
    while len(face_queue) > 0:
        f = face_queue.pop(0)
        halfedges = get_candidate_halfedges(f)
        # Check if halfedge of face is not defined
        # This should only occur if it is the first triangle
        if f.halfedge == None:
            # Should only occur for first triangle
            h = halfedges[0]
            f.halfedge = h
            #h.face = f
        # Set prev and next pointers for all halfedges of face
        h = f.halfedge
        h_start = h
        while True:
            # Set halfedge face
            h.face = f
            # Set up faces of opposite halfedge
            if h.e_op.face == None:
                # Set up op face
                fs = get_candidate_faces(h)
                # If opposite half edge not on boundary
                if len(fs) == 2:
                    other_f = fs[0] if fs[0] != f else fs[1]
                    # Assign face of halfedge
                    h.e_op.face = other_f
                    if other_f.halfedge == None:
                        # Assign halfedge of face
                        other_f.halfedge = h.e_op
                    # Add other_f to queue if in face_set and remove
                    if other_f.boundary == True:
                        face_queue.append(other_f)
                    if other_f in face_set:
                        face_queue.append(other_f)
                        face_set.remove(other_f)
                # Handle case where other_f is boundary.
                elif len(fs) == 1:
                    boundary_halfedges_set.add(h.e_op)
            # Set up prev of halfedge and next of halfedge prev          
            p = get_prev(h, halfedges)
            h.e_prev, p.e_next = p, h
            h = h.e_prev
            # Break loop if we return to h_start
            if h == h_start:
                break
        
    # Handle boundaries
    # ASSUMING NO LONE EDGES
    while len(boundary_halfedges_set) > 0:
        # Get a boundary halfedge.
        h = boundary_halfedges_set.pop()
        # Make a boundary face.
        f = Face(boundary = True)
        # Add boundary face to boundaries list.
        triangulation.boundaries.append(f)
        h_start = h
        h.face = f
        while True:
            # Current boundary edge
            h_cur = h
            # Find next halfedge
            while h.face != None:
                if h.face == f:
                    h = h.e_op.e_prev.e_op
                    if h.face == f:
                        break
                else:
                    h = h.e_prev.e_op
                    if h.face == f:
                        break
            boundary_halfedges_set.discard(h)
            h_cur.e_next, h.e_prev = h, h_cur
            h.face = f
            if h == h_start:
                break
                
    return triangulation
