# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import mock
import numpy as np
import pytest

from tests import utils


# pylint: disable=invalid-name,no-member
slow = pytest.mark.skipif(
    pytest.config.getoption('--ignore-slow') and utils.WITHOUT_SPEEDUPS,
    reason='--ignore-slow ignores the slow tests',
)
# pylint: enable=invalid-name,no-member


class TestSurface(utils.NumPyTestCase):

    REF_TRIANGLE = utils.ref_triangle_uniform_nodes(5)
    REF_TRIANGLE3 = utils.ref_triangle_uniform_nodes(3)
    QUADRATIC = np.asfortranarray([
        [0.0, 0.0],
        [1.25, 0.5],
        [2.0, 1.0],
        [-1.5, 0.75],
        [0.0, 2.0],
        [-3.0, 3.0],
    ])
    UNIT_TRIANGLE = np.asfortranarray([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ])
    ZEROS = np.zeros((3, 2), order='F')

    @staticmethod
    def _get_target_class():
        from bezier import surface

        return surface.Surface

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _make_one_no_slots(self, *args, **kwargs):
        class NoSlots(self._get_target_class()):
            pass

        return NoSlots(*args, **kwargs)

    def test_constructor(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.625, 0.5],
            [1.0, 0.75],
        ])
        surface = self._make_one(nodes, 1, _copy=False)
        self.assertEqual(surface._degree, 1)
        self.assertEqual(surface._dimension, 2)
        self.assertIs(surface._nodes, nodes)
        self.assertEqual(surface._base_x, 0.0)
        self.assertEqual(surface._base_y, 0.0)
        self.assertEqual(surface._width, 1.0)
        self.assertIsNone(surface._area)
        self.assertIsNone(surface._edges)
        self.assertIsNone(surface._is_valid)

    def test_constructor_wrong_dimension(self):
        nodes = np.asfortranarray([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes, 0)

        nodes = np.zeros((3, 2, 2), order='F')
        with self.assertRaises(ValueError):
            self._make_one(nodes, 1)

    def test_from_nodes_factory(self):
        nodes = np.asfortranarray([
            [0.0, 0.0, 0.0],
            [1.0, 0.5, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 2.0],
            [1.0, 1.0, 0.0],
            [2.0, 3.0, 0.0],
        ])
        klass = self._get_target_class()

        surface = klass.from_nodes(
            nodes, base_x=0.25, base_y=0.0, width=0.625)
        self.assertIsInstance(surface, klass)
        self.assertEqual(surface._degree, 2)
        self.assertEqual(surface._dimension, 3)
        self.assertEqual(surface._nodes, nodes)
        self.assertEqual(surface._base_x, 0.25)
        self.assertEqual(surface._base_y, 0.0)
        self.assertEqual(surface._width, 0.625)
        self.assertIsNone(surface._area)
        self.assertIsNone(surface._edges)
        self.assertIsNone(surface._is_valid)

    def test___repr__(self):
        nodes = np.zeros((15, 3), order='F')
        surface = self._make_one(nodes, 4)
        expected = '<Surface (degree=4, dimension=3)>'
        self.assertEqual(repr(surface), expected)

    def test___repr__custom_triangle(self):
        from bezier import surface as surface_mod

        degree = 4
        dimension = 3
        num_nodes = ((degree + 1) * (degree + 2)) // 2
        nodes = np.zeros((num_nodes, dimension), order='F')
        base_x = 0.46875
        base_y = 0.3125
        width = 0.03125
        surface = self._make_one(
            nodes, degree, base_x=base_x, base_y=base_y, width=width)
        expected = surface_mod._REPR_TEMPLATE.format(
            'Surface', degree, dimension, base_x, base_y, width)
        self.assertEqual(repr(surface), expected)

    def test__get_degree_valid(self):
        klass = self._get_target_class()

        self.assertEqual(0, klass._get_degree(1))
        self.assertEqual(1, klass._get_degree(3))
        self.assertEqual(2, klass._get_degree(6))
        self.assertEqual(3, klass._get_degree(10))
        self.assertEqual(11, klass._get_degree(78))

    def test__get_degree_invalid(self):
        klass = self._get_target_class()

        with self.assertRaises(ValueError):
            klass._get_degree(2)

        with self.assertRaises(ValueError):
            klass._get_degree(9)

    def test_area_property_not_cached(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 3.0],
        ])
        surface = self._make_one(nodes, 1)
        self.assertIsNone(surface._area)
        with self.assertRaises(NotImplementedError):
            getattr(surface, 'area')

    def test_area_property(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 3.0],
        ])
        surface = self._make_one(nodes, 1)
        area = 3.14159
        surface._area = area
        self.assertEqual(surface.area, area)

    def test_width_property(self):
        surface = self._make_one(self.ZEROS, 1)
        self.assertEqual(surface.width, 1.0)

    def test_base_x_property(self):
        surface = self._make_one(self.ZEROS, 1)
        self.assertEqual(surface.base_x, 0.0)

    def test_base_y_property(self):
        surface = self._make_one(self.ZEROS, 1)
        self.assertEqual(surface.base_y, 0.0)

    def test__compute_edge_nodes(self):
        nodes = np.asfortranarray([
            [1.0, 2.0],
            [4.0, 2.5],
            [0.0, 4.0],
        ])
        p100, p010, p001 = nodes
        surface = self._make_one(nodes, 1)

        nodes1, nodes2, nodes3 = surface._compute_edge_nodes()
        expected1 = np.asfortranarray(np.vstack([p100, p010]))
        self.assertEqual(nodes1, expected1)
        expected2 = np.asfortranarray(np.vstack([p010, p001]))
        self.assertEqual(nodes2, expected2)
        expected3 = np.asfortranarray(np.vstack([p001, p100]))
        self.assertEqual(nodes3, expected3)

    def _edges_helper(self, edge1, edge2, edge3,
                      nodes1, nodes2, nodes3):
        import bezier

        self.assertIsInstance(edge1, bezier.Curve)
        self.assertEqual(edge1._edge_index, 0)
        self.assertIs(edge1.next_edge, edge2)
        self.assertIs(edge1.previous_edge, edge3)
        self.assertEqual(edge1.nodes, nodes1)

        self.assertIsInstance(edge2, bezier.Curve)
        self.assertEqual(edge2._edge_index, 1)
        self.assertIs(edge2.next_edge, edge3)
        self.assertIs(edge2.previous_edge, edge1)
        self.assertEqual(edge2.nodes, nodes2)

        self.assertIsInstance(edge3, bezier.Curve)
        self.assertEqual(edge3._edge_index, 2)
        self.assertIs(edge3.next_edge, edge1)
        self.assertIs(edge3.previous_edge, edge2)
        self.assertEqual(edge3.nodes, nodes3)

    def test__compute_edges_linear(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [2.0, 1.0],
            [-3.0, 3.0],
        ])
        p100, p010, p001 = nodes
        surface = self._make_one(nodes, 1)

        edge1, edge2, edge3 = surface._compute_edges()
        nodes1 = np.asfortranarray(np.vstack([p100, p010]))
        nodes2 = np.asfortranarray(np.vstack([p010, p001]))
        nodes3 = np.asfortranarray(np.vstack([p001, p100]))
        self._edges_helper(
            edge1, edge2, edge3, nodes1, nodes2, nodes3)

    def test__compute_edges_quadratic(self):
        nodes = self.QUADRATIC
        p200, p110, p020, p101, p011, p002 = nodes
        surface = self._make_one(nodes, 2)

        edge1, edge2, edge3 = surface._compute_edges()
        nodes1 = np.asfortranarray(np.vstack([p200, p110, p020]))
        nodes2 = np.asfortranarray(np.vstack([p020, p011, p002]))
        nodes3 = np.asfortranarray(np.vstack([p002, p101, p200]))
        self._edges_helper(
            edge1, edge2, edge3, nodes1, nodes2, nodes3)

    def test__compute_edges_cubic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.328125, 0.1484375],
            [0.65625, 0.1484375],
            [1.0, 0.0],
            [0.1484375, 0.328125],
            [0.5, 0.5],
            [1.0, 0.53125],
            [0.1484375, 0.65625],
            [0.53125, 1.0],
            [0.0, 1.0],
        ])
        (p300, p210, p120, p030, p201,
         unused_p111, p021, p102, p012, p003) = nodes
        surface = self._make_one(nodes, 3)

        edges = surface._compute_edges()
        self._edges_helper(
            edges[0], edges[1], edges[2],
            np.asfortranarray(np.vstack([p300, p210, p120, p030])),
            np.asfortranarray(np.vstack([p030, p021, p012, p003])),
            np.asfortranarray(np.vstack([p003, p102, p201, p300])))

    def test__get_edges(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)
        compute_mock = mock.Mock(return_value=mock.sentinel.edges)
        surface._compute_edges = compute_mock

        self.assertIsNone(surface._edges)
        self.assertIs(surface._get_edges(), mock.sentinel.edges)
        self.assertIs(surface._edges, mock.sentinel.edges)

        compute_mock.assert_called_once_with()

    def test__get_edges_cached(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)
        compute_mock = mock.Mock()
        surface._compute_edges = compute_mock

        surface._edges = mock.sentinel.edges
        self.assertIs(surface._get_edges(), mock.sentinel.edges)

        compute_mock.assert_not_called()

    def test_edges_property(self):
        nodes = self.UNIT_TRIANGLE
        surface = self._make_one(nodes, 1)

        edge1, edge2, edge3 = surface.edges
        nodes1 = np.asfortranarray(nodes[:2, :])
        nodes2 = np.asfortranarray(nodes[1:, :])
        nodes3 = np.asfortranarray(nodes[(2, 0), :])
        self._edges_helper(edge1, edge2, edge3,
                           nodes1, nodes2, nodes3)

    def test_edges_property_cached(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)

        # Create mock "edges" to be computed.
        sentinel1 = mock.Mock(spec=['_copy'])
        sentinel2 = mock.Mock(spec=['_copy'])
        sentinel3 = mock.Mock(spec=['_copy'])
        expected = sentinel1, sentinel2, sentinel3
        surface._compute_edges = mock.Mock(return_value=expected)

        # Make sure the "edges" when copied just return themselves.
        sentinel1._copy.return_value = sentinel1
        sentinel2._copy.return_value = sentinel2
        sentinel3._copy.return_value = sentinel3

        # Access the property and check the mocks.
        self.assertEqual(surface.edges, expected)
        surface._compute_edges.assert_any_call()
        self.assertEqual(surface._compute_edges.call_count, 1)
        sentinel1._copy.assert_called_once_with()
        sentinel2._copy.assert_called_once_with()
        sentinel3._copy.assert_called_once_with()

        # Access again but make sure no more calls to _compute_edges().
        self.assertEqual(surface.edges, expected)
        self.assertEqual(surface._compute_edges.call_count, 1)

    def test__verify_barycentric(self):
        klass = self._get_target_class()
        # Valid inside.
        self.assertIsNone(klass._verify_barycentric(0.5, 0.25, 0.25))
        # Valid boundary.
        self.assertIsNone(klass._verify_barycentric(0.5, 0.0, 0.5))
        self.assertIsNone(klass._verify_barycentric(0.25, 0.75, 0.0))
        self.assertIsNone(klass._verify_barycentric(0.0, 0.0, 1.0))
        self.assertIsNone(klass._verify_barycentric(0.0, 0.5, 0.5))
        # Invalid sum.
        with self.assertRaises(ValueError):
            klass._verify_barycentric(0.5, 0.5, 0.5)
        # Invalid lamdba1
        with self.assertRaises(ValueError):
            klass._verify_barycentric(-0.5, 0.75, 0.75)
        # Invalid lamdba2.
        with self.assertRaises(ValueError):
            klass._verify_barycentric(0.75, -0.5, 0.75)
        # Invalid lamdba3.
        with self.assertRaises(ValueError):
            klass._verify_barycentric(0.875, 0.25, -0.125)

    def test_evaluate_barycentric(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1, _copy=False)
        lambda_vals = (0.25, 0.0, 0.75)

        # Just make sure we call the helper.
        patch = mock.patch('bezier._surface_helpers.evaluate_barycentric',
                           return_value=mock.sentinel.evaluated)
        with patch as mocked:
            result = surface.evaluate_barycentric(*lambda_vals)
            self.assertIs(result, mock.sentinel.evaluated)
            mocked.assert_called_once_with(
                self.UNIT_TRIANGLE, 1, *lambda_vals)

    def test_evaluate_barycentric_negative_weights_no_verify(self):
        lambda_vals = (0.25, -0.5, 1.25)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.5],
            [0.0, 1.25],
        ])
        surface = self._make_one(nodes, 1)

        self.assertLess(min(lambda_vals), 0.0)
        result = surface.evaluate_barycentric(*lambda_vals, _verify=False)

        expected = np.asfortranarray([[-0.5, 1.3125]])
        self.assertEqual(result, expected)

    def test_evaluate_barycentric_non_unity_weights_no_verify(self):
        lambda_vals = (0.25, 0.25, 0.25)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.5],
            [0.0, 1.25],
        ])
        surface = self._make_one(nodes, 1)

        self.assertNotEqual(sum(lambda_vals), 1.0)
        result = surface.evaluate_barycentric(*lambda_vals, _verify=False)

        expected = np.asfortranarray([[0.25, 0.4375]])
        self.assertEqual(result, expected)

    def test_evaluate_barycentric_multi_wrong_dimension(self):
        surface = self._make_one(self.ZEROS, 1)
        param_vals_1d = np.zeros((4,), order='F')
        with self.assertRaises(ValueError):
            surface.evaluate_barycentric_multi(param_vals_1d)

    def _eval_bary_multi_helper(self, **kwargs):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [2.0, 1.0],
            [-3.0, 2.0],
        ])
        surface = self._make_one(nodes, 1, _copy=False)
        param_vals = np.asfortranarray([[1.0, 0.0, 0.0]])

        patch = mock.patch(
            'bezier._surface_helpers.evaluate_barycentric_multi',
            return_value=mock.sentinel.evaluated)
        with patch as mocked:
            result = surface.evaluate_barycentric_multi(param_vals, **kwargs)
            self.assertEqual(result, mock.sentinel.evaluated)

            mocked.assert_called_once_with(nodes, 1, param_vals, 2)

    def test_evaluate_barycentric_multi(self):
        self._eval_bary_multi_helper()

    def test_evaluate_barycentric_multi_no_verify(self):
        self._eval_bary_multi_helper(_verify=False)

    def test__verify_cartesian(self):
        klass = self._get_target_class()
        # Valid inside.
        self.assertIsNone(klass._verify_cartesian(0.25, 0.25))
        # Valid boundary.
        self.assertIsNone(klass._verify_cartesian(0.0, 0.5))
        self.assertIsNone(klass._verify_cartesian(0.75, 0.0))
        self.assertIsNone(klass._verify_cartesian(0.0, 1.0))
        self.assertIsNone(klass._verify_cartesian(0.5, 0.5))
        # Invalid s.
        with self.assertRaises(ValueError):
            klass._verify_cartesian(-0.5, 0.75)
        # Invalid t.
        with self.assertRaises(ValueError):
            klass._verify_cartesian(0.25, -0.125)
        # Invalid (1 - s - t).
        with self.assertRaises(ValueError):
            klass._verify_cartesian(0.75, 0.75)

    def test_evaluate_cartesian(self):
        s_t_vals = (0.125, 0.125)
        nodes = np.asfortranarray([
            [1.0, 1.0],
            [2.0, 1.5],
            [1.0, 2.75],
        ])
        surface = self._make_one(nodes, 1)

        expected = np.asfortranarray([[1.125, 1.28125]])
        result = surface.evaluate_cartesian(*s_t_vals)
        self.assertEqual(result, expected)

    def test_evaluate_cartesian_no_verify(self):
        s_t_vals = (0.25, 1.0)
        nodes = np.asfortranarray([
            [1.0, 1.0],
            [2.0, 1.5],
            [1.0, 2.75],
        ])
        surface = self._make_one(nodes, 1)

        expected = np.asfortranarray([[1.25, 2.875]])
        result = surface.evaluate_cartesian(*s_t_vals, _verify=False)
        self.assertEqual(result, expected)

    def test_evaluate_cartesian_calls_helper(self):
        nodes = self.ZEROS
        surface = self._make_one_no_slots(nodes, 1, _copy=False)
        patch = mock.patch('bezier._surface_helpers.evaluate_barycentric',
                           return_value=mock.sentinel.point)

        s_val = 0.25
        t_val = 0.25
        with patch as mocked:
            result = surface.evaluate_cartesian(s_val, t_val)
            self.assertIs(result, mock.sentinel.point)

            mocked.assert_called_once_with(nodes, 1, 0.5, s_val, t_val)

    def test_evaluate_cartesian_multi_wrong_dimension(self):
        surface = self._make_one(self.ZEROS, 1)
        param_vals_1d = np.zeros((4,), order='F')
        with self.assertRaises(ValueError):
            surface.evaluate_cartesian_multi(param_vals_1d)

    def _eval_cartesian_multi_helper(self, **kwargs):
        nodes = np.asfortranarray([
            [2.0, 3.0],
            [0.0, 2.0],
            [3.0, 7.5],
        ])
        surface = self._make_one(nodes, 1, _copy=False)
        param_vals = np.asfortranarray([[1.0, 0.0]])

        patch = mock.patch(
            'bezier._surface_helpers.evaluate_cartesian_multi',
            return_value=mock.sentinel.evaluated)
        with patch as mocked:
            result = surface.evaluate_cartesian_multi(param_vals, **kwargs)
            self.assertEqual(result, mock.sentinel.evaluated)

            mocked.assert_called_once_with(nodes, 1, param_vals)

    def test_evaluate_cartesian_multi(self):
        self._eval_cartesian_multi_helper()

    def test_evaluate_cartesian_multi_no_verify(self):
        self._eval_cartesian_multi_helper(_verify=False)

    def test_plot_wrong_dimension(self):
        nodes = np.asfortranarray([
            [0.0, 0.0, 0.0],
            [1.0, 3.0, 4.0],
            [2.0, 6.0, 9.0],
        ])
        surface = self._make_one(nodes, 1, _copy=False)
        with self.assertRaises(NotImplementedError):
            surface.plot(32)

    @mock.patch('bezier._plot_helpers.new_axis')
    @mock.patch('bezier._plot_helpers.add_patch')
    def test_plot_defaults(self, add_patch_mock, new_axis_mock):
        ax = mock.Mock(spec=[])
        new_axis_mock.return_value = ax

        curve = self._make_one(self.UNIT_TRIANGLE, 1, _copy=False)

        pts_per_edge = 16
        result = curve.plot(pts_per_edge)
        self.assertIs(result, ax)

        # Verify mocks.
        new_axis_mock.assert_called_once_with()
        add_patch_mock.assert_called_once_with(
            ax, None, pts_per_edge, *curve._edges)

    @mock.patch('bezier._plot_helpers.new_axis')
    @mock.patch('bezier._plot_helpers.add_patch')
    def test_plot_explicit(self, add_patch_mock, new_axis_mock):
        ax = mock.Mock(spec=['plot'])
        color = (0.5, 0.5, 0.5)
        curve = self._make_one(self.UNIT_TRIANGLE, 1, _copy=False)

        pts_per_edge = 16
        result = curve.plot(
            pts_per_edge, color=color, ax=ax, with_nodes=True)
        self.assertIs(result, ax)

        # Verify mocks.
        new_axis_mock.assert_not_called()
        add_patch_mock.assert_called_once_with(
            ax, color, pts_per_edge, *curve._edges)
        # Check the call to ax.plot(). We can't assert_any_call()
        # since == breaks on NumPy arrays.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        utils.check_plot_call(
            self, call, self.UNIT_TRIANGLE,
            color='black', marker='o', linestyle='None')

    def _subdivide_helper(self, nodes, expected_a, expected_b,
                          expected_c, expected_d):
        klass = self._get_target_class()

        surface = klass.from_nodes(nodes)
        surface_a, surface_b, surface_c, surface_d = surface.subdivide()

        self.assertIsInstance(surface_a, klass)
        self.assertEqual(surface_a._nodes, expected_a)
        self.assertIsInstance(surface_b, klass)
        self.assertEqual(surface_b._nodes, expected_b)
        self.assertIsInstance(surface_c, klass)
        self.assertEqual(surface_c._nodes, expected_c)
        self.assertIsInstance(surface_d, klass)
        self.assertEqual(surface_d._nodes, expected_d)

    def _subdivide_points_check(self, surface):
        # Using the exponent means that we will divide by
        # 2**exp, which can be done without roundoff (for small
        # enough exponents).
        sub_surfaces = surface.subdivide()

        ref_triangle = self.REF_TRIANGLE
        quarter_a = 0.5 * ref_triangle
        quarters = [
            quarter_a,
            np.asfortranarray([0.5, 0.5]) - quarter_a,  # B
            quarter_a + np.asfortranarray([0.5, 0.0]),  # C
            quarter_a + np.asfortranarray([0.0, 0.5]),  # D
        ]

        for sub_surface, quarter in zip(sub_surfaces, quarters):
            # Make sure sub_surface(ref_triangle) == surface(quarter)
            main_vals = surface.evaluate_cartesian_multi(quarter)
            sub_vals = sub_surface.evaluate_cartesian_multi(ref_triangle)
            self.assertEqual(main_vals, sub_vals)

    def test_subdivide_linear(self):
        expected_a = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5],
        ])
        expected_b = np.asfortranarray([
            [0.5, 0.5],
            [0.0, 0.5],
            [0.5, 0.0],
        ])
        expected_c = np.asfortranarray([
            [0.5, 0.0],
            [1.0, 0.0],
            [0.5, 0.5],
        ])
        expected_d = np.asfortranarray([
            [0.0, 0.5],
            [0.5, 0.5],
            [0.0, 1.0],
        ])
        self._subdivide_helper(self.UNIT_TRIANGLE, expected_a,
                               expected_b, expected_c, expected_d)

    @slow
    def test_subdivide_line_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(3, 2), seed=123987, num_bits=8)

        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertEqual(surface.degree, 1)
        self._subdivide_points_check(surface)

    def test_subdivide_quadratic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.25],
            [1.0, 0.0],
            [0.5, 0.75],
            [0.0, 1.0],
            [0.0, 0.5],
        ])
        expected_a = np.asfortranarray([
            [0.0, 0.0],
            [0.25, 0.125],
            [0.5, 0.125],
            [0.25, 0.375],
            [0.25, 0.5],
            [0.25, 0.5],
        ])
        expected_b = np.asfortranarray([
            [0.25, 0.625],
            [0.25, 0.625],
            [0.25, 0.5],
            [0.5, 0.5],
            [0.25, 0.5],
            [0.5, 0.125],
        ])
        expected_c = np.asfortranarray([
            [0.5, 0.125],
            [0.75, 0.125],
            [1.0, 0.0],
            [0.5, 0.5],
            [0.5, 0.5],
            [0.25, 0.625],
        ])
        expected_d = np.asfortranarray([
            [0.25, 0.5],
            [0.25, 0.625],
            [0.25, 0.625],
            [0.25, 0.625],
            [0.0, 0.75],
            [0.0, 0.5],
        ])
        self._subdivide_helper(nodes, expected_a, expected_b,
                               expected_c, expected_d)

    @slow
    def test_subdivide_quadratic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(6, 2), seed=45001, num_bits=8)

        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertEqual(surface.degree, 2)
        self._subdivide_points_check(surface)

    def test_subdivide_cubic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [3.25, 1.5],
            [6.5, 1.5],
            [10.0, 0.0],
            [1.5, 3.25],
            [5.0, 5.0],
            [10.0, 5.25],
            [1.5, 6.5],
            [5.25, 10.0],
            [0.0, 10.0],
        ])
        expected_a = np.asfortranarray([
            [0.0, 0.0],
            [1.625, 0.75],
            [3.25, 1.125],
            [4.90625, 1.125],
            [0.75, 1.625],
            [2.4375, 2.4375],
            [4.3125, 2.875],
            [1.125, 3.25],
            [2.875, 4.3125],
            [1.125, 4.90625],
        ])
        expected_b = np.asfortranarray([
            [6.96875, 6.96875],
            [4.8125, 6.65625],
            [2.875, 5.96875],
            [1.125, 4.90625],
            [6.65625, 4.8125],
            [4.75, 4.75],
            [2.875, 4.3125],
            [5.96875, 2.875],
            [4.3125, 2.875],
            [4.90625, 1.125],
        ])
        expected_c = np.asfortranarray([
            [4.90625, 1.125],
            [6.5625, 1.125],
            [8.25, 0.75],
            [10.0, 0.0],
            [5.96875, 2.875],
            [7.875, 2.9375],
            [10.0, 2.625],
            [6.65625, 4.8125],
            [8.8125, 5.125],
            [6.96875, 6.96875],
        ])
        expected_d = np.asfortranarray([
            [1.125, 4.90625],
            [2.875, 5.96875],
            [4.8125, 6.65625],
            [6.96875, 6.96875],
            [1.125, 6.5625],
            [2.9375, 7.875],
            [5.125, 8.8125],
            [0.75, 8.25],
            [2.625, 10.0],
            [0.0, 10.0],
        ])
        self._subdivide_helper(nodes, expected_a, expected_b,
                               expected_c, expected_d)

    @slow
    def test_subdivide_cubic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(10, 2), seed=346323, num_bits=8)

        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertEqual(surface.degree, 3)
        self._subdivide_points_check(surface)

    @slow
    def test_subdivide_quartic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(15, 2), seed=741002, num_bits=8)

        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertEqual(surface.degree, 4)
        self._subdivide_points_check(surface)

    @slow
    def test_subdivide_on_the_fly(self):
        # Test for a degree where the subdivision is done on the fly
        # rather than via a stored matrix.
        nodes = utils.get_random_nodes(
            shape=(21, 2), seed=446, num_bits=8)
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.

        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertEqual(surface.degree, 5)
        self._subdivide_points_check(surface)

    def test__compute_valid_valid_linear(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1)
        self.assertTrue(surface._compute_valid())

    def test__compute_valid_invalid_linear(self):
        nodes = np.asfortranarray([
            [0.0, 0.0, 0.0],
            [1.0, 2.0, 2.0],
            [2.0, 4.0, 4.0],
        ])
        surface = self._make_one(nodes, 1)
        self.assertFalse(surface._compute_valid())

    def test__compute_valid_quadratic_valid(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, -0.1875],
            [1.0, 0.0],
            [0.1875, 0.5],
            [0.625, 0.625],
            [0.0, 1.0],
        ])
        surface = self._make_one(nodes, 2)
        self.assertTrue(surface._compute_valid())

    def test__compute_valid_quadratic_invalid(self):
        # B(L1, L2, L3) = [L1^2 + L2^2, L2^2 + L3^2]
        nodes = np.asfortranarray([
            [1.0, 0.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
        ])
        surface = self._make_one(nodes, 2)
        self.assertFalse(surface._compute_valid())

    def test__compute_valid_quadratic_bad_dimension(self):
        nodes = np.zeros((6, 3), order='F')
        surface = self._make_one(nodes, 2)
        with self.assertRaises(NotImplementedError):
            surface._compute_valid()

    def test__compute_valid_cubic_valid(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [3.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.25, 1.25],
            [0.0, 2.0],
            [1.25, 2.25],
            [0.0, 3.0],
        ])
        surface = self._make_one(nodes, 3)
        self.assertTrue(surface._compute_valid())

    def test__compute_valid_cubic_invalid(self):
        # B(L1, L2, L3) = [L1^3 + L2^3, L2^3 + L3^3]
        nodes = np.asfortranarray([
            [1.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
        ])
        surface = self._make_one(nodes, 3)
        self.assertFalse(surface._compute_valid())

    def test__compute_valid_bad_degree(self):
        nodes = np.zeros((15, 2), order='F')
        surface = self._make_one(nodes, 4)
        with self.assertRaises(NotImplementedError):
            surface._compute_valid()

    def test_is_valid_property(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1)
        self.assertTrue(surface.is_valid)

    def test_is_valid_property_cached(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)
        compute_valid = mock.Mock()
        surface._compute_valid = compute_valid
        compute_valid.return_value = True

        self.assertIsNone(surface._is_valid)

        # Access the property and check the mocks.
        self.assertTrue(surface.is_valid)
        self.assertTrue(surface._is_valid)
        compute_valid.assert_called_once_with()

        # Access again but make sure no more calls to _compute_valid().
        self.assertTrue(surface.is_valid)
        self.assertEqual(compute_valid.call_count, 1)

    def test___dict___property(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1, _copy=False)
        props_dict = surface.__dict__
        expected = {
            '_nodes': self.UNIT_TRIANGLE,
            '_dimension': 2,
            '_degree': 1,
            '_base_x': 0.0,
            '_base_y': 0.0,
            '_width': 1.0,
            '_area': None,
            '_edges': None,
            '_is_valid': None,
        }
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``surface``.
        expected['_width'] = 1.5
        self.assertNotEqual(surface._width, expected['_width'])

    def test_locate(self):
        surface = self._make_one(self.QUADRATIC, 2)
        point = surface.evaluate_cartesian(0.5, 0.25)
        s, t = surface.locate(point)
        self.assertEqual(s, 0.5)
        self.assertEqual(t, 0.25)

    def test_locate_no_verify(self):
        surface = self._make_one(self.QUADRATIC, 2)
        s = 0.125
        t = 0.125
        x_val, y_val = surface.evaluate_cartesian(s, t).flatten()
        point = np.asfortranarray([
            [x_val, y_val],
            [np.nan, np.nan],
        ])
        # Make sure it fails.
        with self.assertRaises(ValueError):
            surface.locate(point)

        # Will only use the first row if _verify=False.
        computed_s, computed_t = surface.locate(point, _verify=False)
        self.assertEqual(s, computed_s)
        self.assertEqual(t, computed_t)

    def test_locate_bad_dimension(self):
        nodes = np.asfortranarray([[0.0], [1.0], [2.0]])
        surface = self._make_one(nodes, 1)
        with self.assertRaises(NotImplementedError):
            surface.locate(None)

    def test_locate_bad_point(self):
        surface = self._make_one(self.QUADRATIC, 2)
        point1 = np.asfortranarray([0.0, 1.0])
        point2 = np.asfortranarray([[0.0, 1.0, 2.0]])
        with self.assertRaises(ValueError):
            surface.locate(point1)
        with self.assertRaises(ValueError):
            surface.locate(point2)

    def _basic_intersect_helper(self, **kwargs):
        import bezier

        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        # Similar triangle with overlapping square.
        nodes = np.asfortranarray([
            [0.5, 0.0],
            [0.5, 1.0],
            [-0.5, 1.0],
        ])
        surface2 = self._make_one(nodes, 1)

        intersections = surface1.intersect(surface2, **kwargs)
        self.assertEqual(len(intersections), 1)

        intersection = intersections[0]
        self.assertIsInstance(intersection, bezier.CurvedPolygon)

        all_edges = surface1._get_edges() + surface2._get_edges()
        # Check which sides the intersectioned edges come from.
        self.assertEqual(
            [all_edges.index(edge.root) for edge in intersection._edges],
            [5, 3, 1, 2])
        self.assertEqual([edge.start for edge in intersection._edges],
                         [0.5, 0.0, 0.5, 0.0])
        self.assertEqual([edge.end for edge in intersection._edges],
                         [1.0, 0.5, 1.0, 0.5])

    def test_intersect(self):
        self._basic_intersect_helper()

    def test_intersect_no_verify(self):
        self._basic_intersect_helper(_verify=False)

    def test_intersect_algebraic(self):
        from bezier import _intersection_helpers

        strategy = _intersection_helpers.IntersectionStrategy.algebraic
        self._basic_intersect_helper(strategy=strategy)

    def test_intersect_disjoint_bbox(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes = np.asfortranarray([
            [4.0, 0.0],
            [5.0, 0.0],
            [4.0, 1.0],
        ])
        surface2 = self._make_one(nodes, 1)

        intersections = surface1.intersect(surface2)
        self.assertEqual(intersections, [])

    def test_intersect_tangent_bbox(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.0, 1.0],
            [-1.0, 1.0],
        ])
        surface2 = self._make_one(nodes, 1)

        intersections = surface1.intersect(surface2)
        self.assertEqual(intersections, [])

    def test_intersect_non_surface(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1)
        with self.assertRaises(TypeError):
            surface.intersect(object())

    def test_intersect_unsupported_dimension(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes2 = np.zeros((3, 3), order='F')
        surface2 = self._make_one(nodes2, 1)

        with self.assertRaises(NotImplementedError):
            surface1.intersect(surface2)
        with self.assertRaises(NotImplementedError):
            surface2.intersect(surface1)

    def test_elevate_linear(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [2.0, 1.0],
            [-1.0, 2.0],
        ])
        surface = self._make_one(nodes, 1)
        elevated = surface.elevate()
        expected = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.5],
            [2.0, 1.0],
            [-0.5, 1.0],
            [0.5, 1.5],
            [-1.0, 2.0],
        ])
        self.assertEqual(surface.degree, 1)
        self.assertEqual(elevated.degree, 2)
        self.assertEqual(elevated.nodes, expected)

        main_vals = surface.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        sub_vals = elevated.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        self.assertEqual(main_vals, sub_vals)

    def test_elevate_quadratic(self):
        klass = self._get_target_class()

        nodes = np.asfortranarray([[0.0], [6.0], [9.0], [0.0], [6.0], [-3.0]])
        surface = klass.from_nodes(nodes)
        elevated = surface.elevate()
        expected = np.asfortranarray([
            [0.0], [4.0], [7.0], [9.0], [0.0],
            [4.0], [7.0], [-1.0], [3.0], [-3.0]])
        self.assertEqual(surface.degree, 2)
        self.assertEqual(elevated.degree, 3)
        self.assertEqual(elevated.nodes, expected)

        main_vals = surface.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        sub_vals = elevated.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        self.assertEqual(main_vals, sub_vals)
