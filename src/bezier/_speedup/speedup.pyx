#!python
#cython: boundscheck=False, wraparound=False

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


import cython
import numpy as np
cimport numpy as np
from libcpp cimport bool as bool_t


cdef extern from "<stdbool.h>":
    pass

cdef extern:
    void de_casteljau_one_round_fort "de_casteljau_one_round" (
        int num_nodes, int dimension, double *nodes, int degree,
        double lambda1, double lambda2, double lambda3, double *new_nodes)
    void evaluate_curve_barycentric_fort "evaluate_curve_barycentric" (
        int degree, int dimension, double *nodes, int num_vals,
        double *lambda1, double *lambda2, double *evaluated)
    void evaluate_multi_fort "evaluate_multi" (
        int degree, int dimension, double *nodes, int num_vals,
        double *s_vals, double *evaluated)
    void linearization_error_fort "linearization_error" (
        int degree, int dimension, double *nodes, double *error)
    void evaluate_barycentric_fort "evaluate_barycentric" (
        int num_nodes, int dimension, double *nodes, int degree,
        double lambda1, double lambda2, double lambda3, double *point)
    void evaluate_barycentric_multi_fort "evaluate_barycentric_multi" (
        int num_nodes, int dimension, double *nodes, int degree,
        int num_vals, double *param_vals, double *evaluated)
    void evaluate_cartesian_multi_fort "evaluate_cartesian_multi" (
        int num_nodes, int dimension, double *nodes, int degree,
        int num_vals, double *param_vals, double *evaluated)
    void cross_product_fort "cross_product" (double *vec0, double *vec1, double *result)
    void segment_intersection_fort "segment_intersection" (
        double *start0, double *end0, double *start1, double *end1,
        double *s_val, double *s_val, bool_t *success)
    void bbox_fort "bbox" (
        int num_nodes, double *nodes, double *left,
        double *right, double *bottom, double *top)
    void specialize_curve_fort "specialize_curve" (
        int degree, int dimension, double *nodes, double start, double end,
        double curve_start, double curve_end, double *new_nodes,
        double *true_start, double *true_end)
    void jacobian_both_fort "jacobian_both" (
        int num_nodes, int dimension, double *nodes,
        int degree, double *new_nodes)
    void evaluate_hodograph_fort "evaluate_hodograph" (
        double s, int degree, int dimension, double *Nodes, double *hodograph)
    void newton_refine_intersect_fort "newton_refine_intersect" (
        double s, int degree1, double *nodes1,
        double t, int degree2, double *nodes2,
        double *new_s, double *new_t)
    void jacobian_det_fort "jacobian_det" (
        int num_nodes, double *nodes, int degree,
        int num_vals, double *st_vals, double *evaluated)
    void bbox_intersect_fort "bbox_intersect" (
        int num_nodes1, double *nodes1,
        int num_nodes2, double *nodes2, int *enum_val)
    void wiggle_interval_fort "wiggle_interval" (double value, double *result, bool_t *success)
    void parallel_different_fort "parallel_different" (
        double *start0, double *end0,
        double *start1, double *end1, bool_t *result)
    void from_linearized_fort "from_linearized" (
        double error1, double start1, double end1,
        double *start_node1, double *end_node1, int degree1, double *nodes1,
        double error2, double start2, double end2,
        double *start_node2, double *end_node2, int degree2, double *nodes2,
        double *refined_s, double *refined_t,
        bool_t *does_intersect, int *py_exc)


def de_casteljau_one_round(
        np.ndarray[double, ndim=2, mode='fortran'] nodes, int degree,
        double lambda1, double lambda2, double lambda3):
    cdef int num_nodes, dimension
    cdef np.ndarray[double, ndim=2, mode='fortran'] new_nodes

    num_nodes, dimension = np.shape(nodes)
    new_nodes = np.empty((num_nodes - degree - 1, dimension), order='F')

    de_casteljau_one_round_fort(
        num_nodes, dimension, &nodes[0, 0], degree,
        lambda1, lambda2, lambda3, &new_nodes[0, 0])
    return new_nodes


def evaluate_curve_barycentric(
        np.ndarray[double, ndim=2, mode='fortran'] nodes,
        np.ndarray[double, ndim=1, mode='fortran'] lambda1,
        np.ndarray[double, ndim=1, mode='fortran'] lambda2):
    cdef int num_nodes, dimension, degree, num_vals
    cdef np.ndarray[double, ndim=2, mode='fortran'] evaluated

    num_nodes, dimension = np.shape(nodes)
    degree = num_nodes - 1
    num_vals, = np.shape(lambda1)
    evaluated = np.empty((num_vals, dimension), order='F')
    evaluate_curve_barycentric_fort(
        degree,
        dimension,
        &nodes[0, 0],
        num_vals,
        &lambda1[0],
        &lambda2[0],
        &evaluated[0, 0],
    )
    return evaluated


def evaluate_multi(
        np.ndarray[double, ndim=2, mode='fortran'] nodes,
        np.ndarray[double, ndim=1, mode='fortran'] s_vals):
    cdef int num_nodes, dimension, degree, num_vals
    cdef np.ndarray[double, ndim=2, mode='fortran'] evaluated

    num_nodes, dimension = np.shape(nodes)
    degree = num_nodes - 1
    num_vals, = np.shape(s_vals)
    evaluated = np.empty((num_vals, dimension), order='F')
    evaluate_multi_fort(
        degree,
        dimension,
        &nodes[0, 0],
        num_vals,
        &s_vals[0],
        &evaluated[0, 0],
    )
    return evaluated


def linearization_error(np.ndarray[double, ndim=2, mode='fortran'] nodes):
    cdef int num_nodes, dimension, degree
    cdef double error

    num_nodes, dimension = np.shape(nodes)
    degree = num_nodes - 1

    linearization_error_fort(
        degree,
        dimension,
        &nodes[0, 0],
        &error,
    )
    return error


def evaluate_barycentric(
        np.ndarray[double, ndim=2, mode='fortran'] nodes,
        int degree, double lambda1, double lambda2, double lambda3):
    cdef int num_nodes, dimension
    cdef np.ndarray[double, ndim=2, mode='fortran'] point

    num_nodes, dimension = np.shape(nodes)
    point = np.empty((1, dimension), order='F')

    evaluate_barycentric_fort(
        num_nodes,
        dimension,
        &nodes[0, 0],
        degree,
        lambda1,
        lambda2,
        lambda3,
        &point[0, 0],
    )
    return point


def evaluate_barycentric_multi(
        np.ndarray[double, ndim=2, mode='fortran'] nodes, int degree,
        np.ndarray[double, ndim=2, mode='fortran'] param_vals, int dimension):
    cdef int num_nodes, num_vals
    cdef np.ndarray[double, ndim=2, mode='fortran'] evaluated

    # NOTE: We don't check that there are ``dimension`` columns.
    num_nodes, _ = np.shape(nodes)
    # NOTE: We don't check that there are 3 columns.
    num_vals, _ = np.shape(param_vals)
    evaluated = np.empty((num_vals, dimension), order='F')

    evaluate_barycentric_multi_fort(
        num_nodes,
        dimension,
        &nodes[0, 0],
        degree,
        num_vals,
        &param_vals[0, 0],
        &evaluated[0, 0],
    )
    return evaluated


def evaluate_cartesian_multi(
        np.ndarray[double, ndim=2, mode='fortran'] nodes, int degree,
        np.ndarray[double, ndim=2, mode='fortran'] param_vals, int dimension):
    cdef int num_nodes, num_vals
    cdef np.ndarray[double, ndim=2, mode='fortran'] evaluated

    # NOTE: We don't check that there are ``dimension`` columns.
    num_nodes, _ = np.shape(nodes)
    # NOTE: We don't check that there are 2 columns.
    num_vals, _ = np.shape(param_vals)
    evaluated = np.empty((num_vals, dimension), order='F')

    evaluate_cartesian_multi_fort(
        num_nodes,
        dimension,
        &nodes[0, 0],
        degree,
        num_vals,
        &param_vals[0, 0],
        &evaluated[0, 0],
    )
    return evaluated


def cross_product(
        np.ndarray[double, ndim=2, mode='fortran'] vec0,
        np.ndarray[double, ndim=2, mode='fortran'] vec1):
    cdef double result

    cross_product_fort(
        &vec0[0, 0],
        &vec1[0, 0],
        &result,
    )
    return result


def segment_intersection(
        np.ndarray[double, ndim=2, mode='fortran'] start0,
        np.ndarray[double, ndim=2, mode='fortran'] end0,
        np.ndarray[double, ndim=2, mode='fortran'] start1,
        np.ndarray[double, ndim=2, mode='fortran'] end1):
    cdef double s_val, t_val
    cdef bool_t success

    segment_intersection_fort(
        &start0[0, 0],
        &end0[0, 0],
        &start1[0, 0],
        &end1[0, 0],
        &s_val,
        &t_val,
        &success,
    )
    return s_val, t_val, success


def bbox(np.ndarray[double, ndim=2, mode='fortran'] nodes):
    cdef int num_nodes
    cdef double left, right, bottom, top

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bbox_fort(
        num_nodes,
        &nodes[0, 0],
        &left,
        &right,
        &bottom,
        &top,
    )

    return left, right, bottom, top


def specialize_curve(
        np.ndarray[double, ndim=2, mode='fortran'] nodes, double start,
        double end, double curve_start, double curve_end, int degree):
    cdef int num_nodes, dimension
    cdef double true_start, true_end
    cdef np.ndarray[double, ndim=2, mode='fortran'] new_nodes

    num_nodes, dimension = np.shape(nodes)
    new_nodes = np.empty((num_nodes, dimension), order='F')

    specialize_curve_fort(
        degree,
        dimension,
        &nodes[0, 0],
        start,
        end,
        curve_start,
        curve_end,
        &new_nodes[0, 0],
        &true_start,
        &true_end,
    )

    return new_nodes, true_start, true_end


def jacobian_both(
        np.ndarray[double, ndim=2, mode='fortran'] nodes,
        int degree, int dimension):
    cdef int num_nodes
    cdef np.ndarray[double, ndim=2, mode='fortran'] new_nodes

    # NOTE: We don't check that there are ``dimension`` columns.
    num_nodes, _ = np.shape(nodes)
    new_nodes = np.empty((num_nodes - degree - 1, 2 * dimension), order='F')

    jacobian_both_fort(
        num_nodes,
        dimension,
        &nodes[0, 0],
        degree,
        &new_nodes[0, 0],
    )

    return new_nodes


def evaluate_hodograph(
        double s, np.ndarray[double, ndim=2, mode='fortran'] nodes,
        int degree):
    cdef int dimension
    cdef np.ndarray[double, ndim=2, mode='fortran'] hodograph

    # NOTE: We don't check that there are ``degree + 1`` rows.
    _, dimension = np.shape(nodes)

    hodograph = np.empty((1, dimension), order='F')

    evaluate_hodograph_fort(
        s,
        degree,
        dimension,
        &nodes[0, 0],
        &hodograph[0, 0],
    )

    return hodograph


def newton_refine_intersect(
        double s, np.ndarray[double, ndim=2, mode='fortran'] nodes1,
        double t, np.ndarray[double, ndim=2, mode='fortran'] nodes2):
    cdef int num_nodes1, degree1, num_nodes2, degree2
    cdef double new_s, new_t

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    degree1 = num_nodes1 - 1
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(nodes2)
    degree2 = num_nodes2 - 1

    newton_refine_intersect_fort(
        s,
        degree1,
        &nodes1[0, 0],
        t,
        degree2,
        &nodes2[0, 0],
        &new_s,
        &new_t,
    )

    return new_s, new_t


def jacobian_det(
        np.ndarray[double, ndim=2, mode='fortran'] nodes, int degree,
        np.ndarray[double, ndim=2, mode='fortran'] st_vals):
    cdef int num_nodes, num_vals
    cdef np.ndarray[double, ndim=1, mode='fortran'] evaluated

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)
    # NOTE: We don't check that there are 2 columns.
    num_vals, _ = np.shape(st_vals)

    evaluated = np.empty((num_vals,), order='F')

    jacobian_det_fort(
        num_nodes,
        &nodes[0, 0],
        degree,
        num_vals,
        &st_vals[0, 0],
        &evaluated[0],
    )

    return evaluated


def bbox_intersect(
        np.ndarray[double, ndim=2, mode='fortran'] nodes1,
        np.ndarray[double, ndim=2, mode='fortran'] nodes2):
    cdef int num_nodes1, num_nodes2, enum_val

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(nodes2)

    bbox_intersect_fort(
        num_nodes1,
        &nodes1[0, 0],
        num_nodes2,
        &nodes2[0, 0],
        &enum_val,
    )

    return enum_val


def wiggle_interval(double value):
    cdef double result
    cdef bool_t success

    wiggle_interval_fort(
        value,
        &result,
        &success,
    )
    return result, success


def parallel_different(
        np.ndarray[double, ndim=2, mode='fortran'] start0,
        np.ndarray[double, ndim=2, mode='fortran'] end0,
        np.ndarray[double, ndim=2, mode='fortran'] start1,
        np.ndarray[double, ndim=2, mode='fortran'] end1):
    cdef bool_t result

    parallel_different_fort(
        &start0[0, 0],
        &end0[0, 0],
        &start1[0, 0],
        &end1[0, 0],
        &result,
    )
    return result


def from_linearized(
        double error1, double start1, double end1,
        np.ndarray[double, ndim=2, mode='fortran'] start_node1,
        np.ndarray[double, ndim=2, mode='fortran'] end_node1,
        np.ndarray[double, ndim=2, mode='fortran'] nodes1,
        double error2, double start2, double end2,
        np.ndarray[double, ndim=2, mode='fortran'] start_node2,
        np.ndarray[double, ndim=2, mode='fortran'] end_node2,
        np.ndarray[double, ndim=2, mode='fortran'] nodes2):
    cdef int num_nodes1, degree1, num_nodes2, degree2
    cdef double refined_s, refined_t
    cdef bool_t does_intersect
    cdef int py_exc

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    degree1 = num_nodes1 - 1
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(nodes2)
    degree2 = num_nodes2 - 1

    from_linearized_fort(
        error1,
        start1,
        end1,
        &start_node1[0, 0],
        &end_node1[0, 0],
        degree1,
        &nodes1[0, 0],
        error2,
        start2,
        end2,
        &start_node2[0, 0],
        &end_node2[0, 0],
        degree2,
        &nodes2[0, 0],
        &refined_s,
        &refined_t,
        &does_intersect,
        &py_exc,
    )

    if py_exc == 1:
        raise NotImplementedError('Line segments parallel.')
    elif py_exc == 2:
        raise ValueError('outside of unit interval')

    return refined_s, refined_t, does_intersect
