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

"""``ctypes`` interface to ``libspeedup`` shared library."""


import ctypes
import os

import numpy as np


HERE = os.path.abspath(os.path.dirname(__file__))
SHARED_LIB = os.path.join(HERE, 'libspeedup.so')
if os.path.exists(SHARED_LIB):
    LIBSPEEDUP = ctypes.cdll.LoadLibrary(SHARED_LIB)
else:
    LIBSPEEDUP = None
    MSG = '{!r} is missing'.format(SHARED_LIB)
    raise ImportError(MSG)


def numpy_pointer(array):
    assert array.flags.f_contiguous
    return array.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def de_casteljau_one_round(nodes, degree, lambda1, lambda2, lambda3):
    num_nodes, dimension = nodes.shape
    new_nodes = np.empty((num_nodes - degree - 1, dimension), order='F')
    LIBSPEEDUP.de_casteljau_one_round(
        ctypes.c_int(num_nodes),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(degree),
        ctypes.c_double(lambda1),
        ctypes.c_double(lambda2),
        ctypes.c_double(lambda3),
        numpy_pointer(new_nodes),
    )
    return new_nodes


def evaluate_curve_barycentric(nodes, lambda1, lambda2):
    num_nodes, dimension = nodes.shape
    degree = num_nodes - 1
    num_vals, = lambda1.shape
    evaluated = np.empty((num_vals, dimension), order='F')
    LIBSPEEDUP.evaluate_curve_barycentric(
        ctypes.c_int(degree),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(num_vals),
        numpy_pointer(lambda1),
        numpy_pointer(lambda2),
        numpy_pointer(evaluated),
    )
    return evaluated


def evaluate_multi(nodes, s_vals):
    num_nodes, dimension = nodes.shape
    degree = num_nodes - 1
    num_vals, = s_vals.shape
    evaluated = np.empty((num_vals, dimension), order='F')
    LIBSPEEDUP.evaluate_multi(
        ctypes.c_int(degree),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(num_vals),
        numpy_pointer(s_vals),
        numpy_pointer(evaluated),
    )
    return evaluated


def linearization_error(nodes):
    num_nodes, dimension = nodes.shape
    degree = num_nodes - 1
    error = ctypes.c_double()

    LIBSPEEDUP.linearization_error(
        ctypes.c_int(degree),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.pointer(error),
    )
    return error.value


def evaluate_barycentric(nodes, degree, lambda1, lambda2, lambda3):
    num_nodes, dimension = nodes.shape
    point = np.empty((1, dimension), order='F')

    LIBSPEEDUP.evaluate_barycentric(
        ctypes.c_int(num_nodes),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(degree),
        ctypes.c_double(lambda1),
        ctypes.c_double(lambda2),
        ctypes.c_double(lambda3),
        numpy_pointer(point),
    )
    return point


def evaluate_barycentric_multi(nodes, degree, param_vals):
    num_nodes, dimension = nodes.shape
    num_vals, cols = param_vals.shape
    assert cols == 3
    evaluated = np.empty((num_vals, dimension), order='F')

    LIBSPEEDUP.evaluate_barycentric_multi(
        ctypes.c_int(num_nodes),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(degree),
        ctypes.c_int(num_vals),
        numpy_pointer(param_vals),
        numpy_pointer(evaluated),
    )
    return evaluated


def evaluate_cartesian_multi(nodes, degree, param_vals):
    num_nodes, dimension = nodes.shape
    num_vals, cols = param_vals.shape
    assert cols == 2
    evaluated = np.empty((num_vals, dimension), order='F')

    LIBSPEEDUP.evaluate_cartesian_multi(
        ctypes.c_int(num_nodes),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(degree),
        ctypes.c_int(num_vals),
        numpy_pointer(param_vals),
        numpy_pointer(evaluated),
    )
    return evaluated


def cross_product(vec0, vec1):
    result = ctypes.c_double()
    LIBSPEEDUP.cross_product(
        numpy_pointer(vec0),
        numpy_pointer(vec1),
        ctypes.pointer(result),
    )
    return result.value


def segment_intersection(start0, end0, start1, end1):
    s_val = ctypes.c_double()
    t_val = ctypes.c_double()
    success = ctypes.c_bool()
    LIBSPEEDUP.segment_intersection(
        numpy_pointer(start0),
        numpy_pointer(end0),
        numpy_pointer(start1),
        numpy_pointer(end1),
        ctypes.pointer(s_val),
        ctypes.pointer(t_val),
        ctypes.pointer(success),
    )
    return s_val.value, t_val.value, success.value


def bbox(nodes):
    num_nodes, cols = nodes.shape
    assert cols == 2
    left = ctypes.c_double()
    right = ctypes.c_double()
    bottom = ctypes.c_double()
    top = ctypes.c_double()

    LIBSPEEDUP.bbox(
        ctypes.c_int(num_nodes),
        numpy_pointer(nodes),
        ctypes.pointer(left),
        ctypes.pointer(right),
        ctypes.pointer(bottom),
        ctypes.pointer(top),
    )

    return left.value, right.value, bottom.value, top.value


def specialize_curve(nodes, start, end, curve_start, curve_end):
    num_nodes, dimension = nodes.shape
    degree = num_nodes - 1
    new_nodes = np.empty((num_nodes, dimension), order='F')

    true_start = ctypes.c_double()
    true_end = ctypes.c_double()

    LIBSPEEDUP.specialize_curve(
        ctypes.c_int(degree),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_double(start),
        ctypes.c_double(end),
        ctypes.c_double(curve_start),
        ctypes.c_double(curve_end),
        numpy_pointer(new_nodes),
        ctypes.pointer(true_start),
        ctypes.pointer(true_end),
    )

    return new_nodes, true_start.value, true_end.value


def jacobian_both(nodes, degree, dimension):
    num_nodes, dimension = nodes.shape
    new_nodes = np.empty((num_nodes - degree - 1, 2 * dimension), order='F')

    LIBSPEEDUP.jacobian_both(
        ctypes.c_int(num_nodes),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        ctypes.c_int(degree),
        numpy_pointer(new_nodes),
    )

    return new_nodes


def evaluate_hodograph(s, nodes):
    num_nodes, dimension = nodes.shape
    degree = num_nodes - 1
    hodograph = np.empty((1, dimension), order='F')

    LIBSPEEDUP.evaluate_hodograph(
        ctypes.c_double(s),
        ctypes.c_int(degree),
        ctypes.c_int(dimension),
        numpy_pointer(nodes),
        numpy_pointer(hodograph),
    )

    return hodograph


def newton_refine_intersect(s, nodes1, t, nodes2):
    num_nodes1, cols = nodes1.shape
    assert cols == 2
    degree1 = num_nodes1 - 1

    num_nodes2, cols = nodes2.shape
    assert cols == 2
    degree2 = num_nodes2 - 1

    new_s = ctypes.c_double()
    new_t = ctypes.c_double()

    LIBSPEEDUP.newton_refine_intersect(
        ctypes.c_double(s),
        ctypes.c_int(degree1),
        numpy_pointer(nodes1),
        ctypes.c_double(t),
        ctypes.c_int(degree2),
        numpy_pointer(nodes2),
        ctypes.pointer(new_s),
        ctypes.pointer(new_t),
    )

    return new_s.value, new_t.value


def jacobian_det(nodes, degree, st_vals):
    num_nodes, cols = nodes.shape
    assert cols == 2

    num_vals, cols = st_vals.shape
    assert cols == 2

    evaluated = np.empty((num_vals,), order='F')

    LIBSPEEDUP.jacobian_det(
        ctypes.c_int(num_nodes),
        numpy_pointer(nodes),
        ctypes.c_int(degree),
        ctypes.c_int(num_vals),
        numpy_pointer(st_vals),
        numpy_pointer(evaluated),
    )

    return evaluated


def bbox_intersect(nodes1, nodes2):
    num_nodes1, cols = nodes1.shape
    assert cols == 2
    num_nodes2, cols = nodes2.shape
    assert cols == 2

    enum_val = ctypes.c_int()

    LIBSPEEDUP.bbox_intersect(
        ctypes.c_int(num_nodes1),
        numpy_pointer(nodes1),
        ctypes.c_int(num_nodes2),
        numpy_pointer(nodes2),
        ctypes.pointer(enum_val),
    )

    return enum_val.value


def wiggle_interval(value):
    result = ctypes.c_double()
    success = ctypes.c_bool()
    LIBSPEEDUP.wiggle_interval(
        ctypes.c_double(value),
        ctypes.pointer(result),
        ctypes.pointer(success),
    )
    return result.value, success.value


def parallel_different(start0, end0, start1, end1):
    result = ctypes.c_bool()
    LIBSPEEDUP.parallel_different(
        numpy_pointer(start0),
        numpy_pointer(end0),
        numpy_pointer(start1),
        numpy_pointer(end1),
        ctypes.pointer(result),
    )
    return result.value


def from_linearized(
        error1, start1, end1, start_node1, end_node1, nodes1,
        error2, start2, end2, start_node2, end_node2, nodes2):
    num_nodes1, cols = nodes1.shape
    assert cols == 2
    degree1 = num_nodes1 - 1

    num_nodes2, cols = nodes2.shape
    assert cols == 2
    degree2 = num_nodes2 - 1

    refined_s = ctypes.c_double()
    refined_t = ctypes.c_double()
    does_intersect = ctypes.c_bool()
    py_exc = ctypes.c_int()

    LIBSPEEDUP.from_linearized(
        ctypes.c_double(error1),
        ctypes.c_double(start1),
        ctypes.c_double(end1),
        numpy_pointer(start_node1),
        numpy_pointer(end_node1),
        ctypes.c_int(degree1),
        numpy_pointer(nodes1),
        ctypes.c_double(error2),
        ctypes.c_double(start2),
        ctypes.c_double(end2),
        numpy_pointer(start_node2),
        numpy_pointer(end_node2),
        ctypes.c_int(degree2),
        numpy_pointer(nodes2),
        ctypes.pointer(refined_s),
        ctypes.pointer(refined_t),
        ctypes.pointer(does_intersect),
        ctypes.pointer(py_exc),
    )

    if py_exc.value == 1:
        raise NotImplementedError('Line segments parallel.')
    elif py_exc.value == 2:
        raise ValueError('outside of unit interval')

    return refined_s.value, refined_t.value, does_intersect.value
