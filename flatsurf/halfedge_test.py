import numpy as np
import numpy.testing as npt
import unittest
import flatsurf.halfedge as halfedge


class VertexTest(unittest.TestCase):
    def test_coord_2d(self):
        coord_2d_orig = np.random.rand(3, 1)
        coord_2d_prime = np.random.rand(3, 1)
        vertex_0 = halfedge.Vertex(coord_2d_orig=coord_2d_orig,
                                   coord_2d_prime=None)
        vertex_1 = halfedge.Vertex(coord_2d_orig=coord_2d_orig,
                                   coord_2d_prime=coord_2d_prime)
        npt.assert_allclose(vertex_0.coord_2d, coord_2d_orig, atol=1e-5)
        npt.assert_allclose(vertex_1.coord_2d, coord_2d_prime, atol=1e-5)

    def test_get_halfedges(self):
        pass

    def test_get_faces(self):
        pass


class FaceTest(unittest.TestCase):
    def test_get_vertices(self):
        pass

    def test_get_halfedges(self):
        pass

    def test_get_adjacent_faces(self):
        pass

    def test_get_shared_faces(self):
        pass

    def test_get_shared_edge(self):
        pass


class HalfEdgeTriangulationTest(unittest.TestCase):
    def test_create_edge_set_and_faces_from_simplices(self):

        coords = np.random.randn(3,3).tolist()
        simplices = [[0, 1, 2]]

        triangulation = halfedge.HalfEdgeTriangulation([], [], None)
        triangulation.vertices = [halfedge.Vertex(coord_3d=coord) for coord in coords]

        # Initialize connectivity dictionaries.
        triangulation._create_connectivity_dicts()
        # Create edge set and faces from simplices
        triangulation._create_edge_set_and_faces_from_simplices(simplices)

        self.assertEqual(len(triangulation.vertices), 3)
        self.assertEqual(len(triangulation.edge_set), 3)
        self.assertEqual(len(triangulation.faces), 1)

        face = triangulation.faces[0]
        for vertex in triangulation.vertices:
            self.assertEqual(triangulation.vertex_to_face_dict[vertex]['faces'][0], face)


    def test_create_halfedges(self):

        coords = np.random.randn(4, 3).tolist()
        simplices = [[0, 1, 2], [1, 2, 3]]

        triangulation = halfedge.HalfEdgeTriangulation([], [], None)
        triangulation.vertices = [halfedge.Vertex(coord_3d=coord) for coord in coords]

        # Initialize connectivity dictionaries.
        triangulation._create_connectivity_dicts()
        # Create edge set and faces from simplices
        triangulation._create_edge_set_and_faces_from_simplices(simplices)

        triangulation._create_halfedges()

        self.assertEqual(len(
            triangulation.vertex_to_face_dict[triangulation.vertices[0]]['halfedges']),
            2
        )

        self.assertEqual(len(
            triangulation.vertex_to_face_dict[triangulation.vertices[1]]['halfedges']),
            3
        )

        self.assertEqual(len(
            triangulation.vertex_to_face_dict[triangulation.vertices[2]]['halfedges']),
            3
        )

        self.assertEqual(len(
            triangulation.vertex_to_face_dict[triangulation.vertices[3]]['halfedges']),
            2
        )

    def test_get_face_candidate_halfedges(self):
        coords = np.random.randn(4, 3).tolist()
        simplices = [[0, 1, 2], [1, 2, 3]]

        triangulation = halfedge.HalfEdgeTriangulation([], [], None)
        triangulation.vertices = [halfedge.Vertex(coord_3d=coord) for coord in coords]

        # Initialize connectivity dictionaries.
        triangulation._create_connectivity_dicts()
        # Create edge set and faces from simplices
        triangulation._create_edge_set_and_faces_from_simplices(simplices)

        triangulation._create_halfedges()

        for face in triangulation.faces:
            # Check that ValueError isn't thrown.
            triangulation._get_face_candidate_halfedges(face)

    def test_create_triangulation_from_coords_and_simplices(self):
        # One Triangle
        a = halfedge.HalfEdgeTriangulation.create_triangulation_from_coords_and_simplices(
                [np.array([0, 0, 0]), np.array([1, 1, 1]),
                 np.array([1, 0, 1])], [[0, 1, 2]])
        self.assertEqual(len(a.boundaries), 1)
        self.assertEqual(len(a.faces), 1)

        # Two triangle
        b = halfedge.HalfEdgeTriangulation.create_triangulation_from_coords_and_simplices(
            [np.array([0, 0, 0]), np.array([1, 1, 1]), np.array([1, 0, 1]), np.array([1, 0, 0])],
            [[0, 1, 2], [0, 1, 3]])
        self.assertEqual(len(b.boundaries), 1)
        self.assertEqual(len(b.faces), 2)

        # Two boundaries
        c = halfedge.HalfEdgeTriangulation.create_triangulation_from_coords_and_simplices(
            [np.array([0, 1, 1]), np.array([0, -1, 1]), np.array([0, -1, -1]), np.array([0, 1, -1]),
             np.array([0, 2, 2]), np.array([0, -2, 2]), np.array([0, -2, -2]), np.array([0, 2, -2])],
            [[0, 1, 4], [1, 4, 5], [1, 2, 5], [2, 5, 6], [2, 3, 6], [3, 6, 7], [0, 3, 7], [0, 4, 7]]
        )
        self.assertEqual(len(c.boundaries), 2)
        self.assertEqual(len(c.faces), 8)

if __name__ == '__main__':
    unittest.main()