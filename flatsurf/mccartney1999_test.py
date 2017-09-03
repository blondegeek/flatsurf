import numpy as np
import numpy.testing as npt
import unittest
import flatsurf.halfedge as halfedge
import flatsurf.mccartney1999 as mc1999

class McCartnet1999FlatteningTest(unittest.TestCase):

    def setUp(self):
        coords = [np.array([0., 0., 0.]), np.array([1., 1., 0.]),
                  np.array([1., 0., 0.]), np.array([0., 1., 0.])]
        simplices = [[0, 1, 2], [0, 1, 3]]
        self.triangulation = halfedge.HalfEdgeTriangulation.from_coords_and_simplices(coords, simplices)

    def test_get_seed_face(self):
        flat_trig = mc1999.McCartney1999Flattening(self.triangulation)
        seed_face = flat_trig.get_seed_face()
        self.assertEqual(type(seed_face), halfedge.Face)

    def test_flatten(self):
        flat_trig = mc1999.McCartney1999Flattening(self.triangulation)
        flat_trig.flatten()
        vertices = flat_trig.vertices
        # Check that trivial triangulation is flattened properly.
        npt.assert_allclose(np.linalg.norm(vertices[0].coord_2d - vertices[1].coord_2d),
                            np.sqrt(2))
        npt.assert_allclose(np.linalg.norm(vertices[2].coord_2d - vertices[3].coord_2d),
                            np.sqrt(2))
        npt.assert_allclose(np.linalg.norm(vertices[0].coord_2d - vertices[2].coord_2d),
                            1.)
        npt.assert_allclose(np.linalg.norm(vertices[0].coord_2d - vertices[3].coord_2d),
                            1.)

class GetCircleIntersectionTest(unittest.TestCase):
    def test_get_circles_intersection(self):
        center1 = np.array([0., 0.])
        center2 = np.array([0., 1.])
        radius1 = 0.5
        radius2 = 0.5
        npt.assert_allclose(mc1999.get_circles_intersection(center1, radius1, center2, radius2)[0],
                            np.array([0, 0.5]))

if __name__ == '__main__':
    unittest.main()