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

"""Private helper methods for B |eacute| zier curves.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import functools

import numpy as np
import six

try:
    import scipy.integrate as _scipy_int
except ImportError:  # pragma: NO COVER
    _scipy_int = None

from bezier import _helpers
try:
    from bezier import _speedup
    import bezier._speedup.speedup  # noqa: F401
except ImportError:  # pragma: NO COVER
    _speedup = None


_MAX_LOCATE_SUBDIVISIONS = 20
_LOCATE_STD_CAP = 0.5**20
_FLOAT64 = np.float64  # pylint: disable=no-member
_REDUCE_THRESHOLD = 0.5**26  # sqrt(machine precision)
# Projections onto the space of degree-elevated nodes.
# If v --> Ev is the elevation map, then P = E (E^T E)^{-1} E^T
# is the projection.
# pylint: disable=bad-whitespace
_PROJECTION0 = np.asfortranarray([
    [0.5, 0.5],
    [0.5, 0.5],
])
_PROJ_DENOM0 = 1.0
_PROJECTION1 = np.asfortranarray([
    [ 2.5, 1.0, -0.5],  # noqa: E201
    [ 1.0, 1.0,  1.0],  # noqa: E201
    [-0.5, 1.0,  2.5],
])
_PROJ_DENOM1 = 3.0
_PROJECTION2 = np.asfortranarray([
    [ 4.75,  0.75, -0.75,  0.25],  # noqa: E201
    [ 0.75,  2.75,  2.25, -0.75],  # noqa: E201
    [-0.75,  2.25,  2.75,  0.75],
    [ 0.25, -0.75,  0.75,  4.75],  # noqa: E201
])
_PROJ_DENOM2 = 5.0
_PROJECTION3 = np.asfortranarray([
    [103.5,   6.0, -9.0,   6.0,  -1.5],
    [  6.0,  81.0, 36.0, -24.0,   6.0],  # noqa: E201
    [ -9.0,  36.0, 51.0,  36.0,  -9.0],  # noqa: E201
    [  6.0, -24.0, 36.0,  81.0,   6.0],  # noqa: E201
    [ -1.5,   6.0, -9.0,   6.0, 103.5],  # noqa: E201
])
_PROJ_DENOM3 = 105.0
# Reductions for a set of degree-elevated nodes.
# If v --> Ev is the elevation map, then R = (E^T E)^{-1} E^T -- the
# pseudo-inverse of E -- actually reduces a set of nodes.
_REDUCTION0 = np.asfortranarray([
    [0.5, 0.5],
])
_REDUCTION_DENOM0 = 1.0
_REDUCTION1 = np.asfortranarray([
    [ 2.5, 1.0, -0.5],  # noqa: E201
    [-0.5, 1.0,  2.5],
])
_REDUCTION_DENOM1 = 3.0
_REDUCTION2 = np.asfortranarray([
    [ 4.75,  0.75, -0.75,  0.25],  # noqa: E201
    [-1.25,  3.75,  3.75, -1.25],
    [ 0.25, -0.75,  0.75,  4.75],  # noqa: E201
])
_REDUCTION_DENOM2 = 5.0
_REDUCTION3 = np.asfortranarray([
    [103.5,   6.0, -9.0,   6.0,  -1.5],
    [-26.5, 106.0, 51.0, -34.0,   8.5],
    [  8.5, -34.0, 51.0, 106.0, -26.5],  # noqa: E201
    [ -1.5,   6.0, -9.0,   6.0, 103.5],  # noqa: E201
])
_REDUCTION_DENOM3 = 105.0
# pylint: enable=bad-whitespace


def make_subdivision_matrices(degree):
    """Make the matrix used to subdivide a curve.

    Args:
        degree (int): The degree of the curve.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray]: The matrices used to convert
           the nodes into left and right nodes, respectively.
    """
    left = np.zeros((degree + 1, degree + 1), order='F')
    right = np.zeros((degree + 1, degree + 1), order='F')
    left[0, 0] = 1.0
    right[-1, -1] = 1.0
    for row in six.moves.xrange(1, degree + 1):
        half_prev = 0.5 * left[row - 1, :row]
        left[row, :row] = half_prev
        left[row, 1:row + 1] += half_prev
        # Populate the complement row (in right) as well.
        complement = degree - row
        # NOTE: We "should" reverse the results when using
        #       the complement, but they are symmetric so
        #       that would be a waste.
        right[complement, -(row + 1):] = left[row, :row + 1]

    return left, right


def _evaluate_multi(nodes, s_vals):
    r"""Computes multiple points along a curve.

    Does so via a modified Horner's method for each value in ``s_vals``
    rather than using the de Casteljau algorithm.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        s_vals (numpy.ndarray): Parameters along the curve (as a
            1D array).

    Returns:
        numpy.ndarray: The evaluated points on the curve as a two dimensional
        NumPy array, with the rows corresponding to each ``s``
        value and the columns to the dimension.
    """
    one_less = 1.0 - s_vals
    return evaluate_multi_barycentric(nodes, one_less, s_vals)


def _evaluate_multi_barycentric(nodes, lambda1, lambda2):
    r"""Evaluates a B |eacute| zier type-function.

    Of the form

    .. math::

       B(\lambda_1, \lambda_2) = \sum_j \binom{n}{j}
           \lambda_1^{n - j} \lambda_2^j \cdot v_j

    for some set of vectors :math:`v_j` given by ``nodes``.

    Does so via a modified Horner's method for each pair of values
    in ``lambda1`` and ``lambda2``, rather than using the
    de Casteljau algorithm.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        lambda1 (numpy.ndarray): Parameters along the curve (as a
            1D array).
        lambda2 (numpy.ndarray): Parameters along the curve (as a
            1D array). Typically we have ``lambda1 + lambda2 == 1``.

    Returns:
        numpy.ndarray: The evaluated points as a two dimensional
        NumPy array, with the rows corresponding to each pair of parameter
        values and the columns to the dimension.
    """
    # NOTE: We assume but don't check that lambda2 has the same shape.
    num_vals, = lambda1.shape
    num_nodes, dimension = nodes.shape
    degree = num_nodes - 1

    # Resize for broadcasting.
    lambda1 = lambda1[:, np.newaxis]
    lambda2 = lambda2[:, np.newaxis]

    result = np.zeros((num_vals, dimension), order='F')
    result += lambda1 * nodes[0, :]

    binom_val = 1.0
    lambda2_pow = np.ones((num_vals, 1), order='F')
    for index in six.moves.xrange(1, degree):
        lambda2_pow *= lambda2
        binom_val = (binom_val * (degree - index + 1)) / index
        result += binom_val * lambda2_pow * nodes[index, :]
        result *= lambda1

    result += lambda2 * lambda2_pow * nodes[degree, :]

    return result


def _vec_size(nodes, s_val):
    r"""Compute :math:`\|B(s)\|_2`.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        s_val (float): Parameter to compute :math:`B(s)`.

    Returns:
        float: The norm of :math:`B(s)`.
    """
    result_vec = evaluate_multi(nodes, np.asfortranarray([s_val]))
    # NOTE: We convert to 1D to make sure NumPy uses vector norm.
    return np.linalg.norm(result_vec[0, :], ord=2)


def compute_length(nodes, degree):
    r"""Approximately compute the length of a curve.

    .. _QUADPACK: https://en.wikipedia.org/wiki/QUADPACK

    If ``degree`` is :math:`n`, then the Hodograph curve
    :math:`B'(s)` is degree :math:`d = n - 1`. Using this curve, we
    approximate the integral:

    .. math::

       \ell\left(B\right) =
           \int_0^1 \| B'(s) \|_2 \, ds

    using `QUADPACK`_ (via SciPy).

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.

    Returns:
        float: The length of the curve.

    Raises:
        OSError: If SciPy is not installed.
    """
    # NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    #       here. This is so we don't re-compute the nodes for the first
    #       derivative every time it is evaluated.
    first_deriv = degree * (nodes[1:, :] - nodes[:-1, :])
    if degree == 1:
        # NOTE: We convert to 1D to make sure NumPy uses vector norm.
        return np.linalg.norm(first_deriv[0, :], ord=2)

    if _scipy_int is None:
        raise OSError('This function requires SciPy for quadrature.')

    size_func = functools.partial(_vec_size, first_deriv)
    length, _ = _scipy_int.quad(size_func, 0.0, 1.0)
    return length


def elevate_nodes(nodes, degree, dimension):
    r"""Degree-elevate a B |eacute| zier curves.

    Does this by converting the current nodes :math:`v_0, \ldots, v_n`
    to new nodes :math:`w_0, \ldots, w_{n + 1}` where

    .. math::

       \begin{align*}
       w_0 &= v_0 \\
       w_j &= \frac{j}{n + 1} v_{j - 1} + \frac{n + 1 - j}{n + 1} v_j \\
       w_{n + 1} &= v_n
       \end{align*}

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.
        dimension (int): The dimension of the curve.

    Returns:
        numpy.ndarray: The nodes of the degree-elevated curve.
    """
    new_nodes = np.zeros((degree + 2, dimension), order='F')

    multipliers = np.arange(1, degree + 1, dtype=_FLOAT64)[:, np.newaxis]
    denominator = degree + 1.0
    new_nodes[1:-1, :] = (
        multipliers * nodes[:-1, :] +
        (denominator - multipliers) * nodes[1:, :])
    # Hold off on division until the end, to (attempt to) avoid round-off.
    new_nodes /= denominator

    # After setting the internal nodes (which require division), set the
    # boundary nodes.
    new_nodes[0, :] = nodes[0, :]
    new_nodes[-1, :] = nodes[-1, :]

    return new_nodes


def de_casteljau_one_round(nodes, lambda1, lambda2):
    """Perform one round of de Casteljau's algorithm.

    The weights are assumed to sum to one.

    Args:
        nodes (numpy.ndarray): Control points for a curve.
        lambda1 (float): First barycentric weight on interval.
        lambda2 (float): Second barycentric weight on interval.

    Returns:
        numpy.ndarray: The nodes for a "blended" curve one degree
        lower.
    """
    return np.asfortranarray(
        lambda1 * nodes[:-1, :] + lambda2 * nodes[1:, :])


def _specialize_curve(nodes, start, end, curve_start, curve_end):
    """Specialize a curve to a re-parameterization

    Does so by taking two points along the number line and then
    reparameterizing the curve onto the interval formed by the
    start and end points.

    .. note::

       This assumes the curve is degree 1 or greater but doesn't check.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Control points for a curve.
        start (float): The start point of the interval we are specializing to.
        end (float): The end point of the interval we are specializing to.
        curve_start (float): The beginning of the original interval which the
            the curve defined by ``nodes`` defines.
        curve_end (float): The end of the original interval which the
            the curve defined by ``nodes`` defines.

    Returns:
        Tuple[numpy.ndarray, float, float]: The control points for the
        specialized curve, as well as the true start and true end of the
        newly created curve.
    """
    # pylint: disable=too-many-locals
    num_nodes, _ = nodes.shape
    degree = num_nodes - 1

    # Uses start-->0, end-->1 to represent the specialization used.
    weights = (
        (1.0 - start, start),
        (1.0 - end, end),
    )
    partial_vals = {
        (0,): de_casteljau_one_round(nodes, *weights[0]),
        (1,): de_casteljau_one_round(nodes, *weights[1]),
    }

    for _ in six.moves.xrange(degree - 1, 0, -1):
        new_partial = {}
        for key, sub_nodes in six.iteritems(partial_vals):
            # Our keys are ascending so we increment from the last value.
            for next_id in six.moves.xrange(key[-1], 1 + 1):
                new_key = key + (next_id,)
                new_partial[new_key] = de_casteljau_one_round(
                    sub_nodes, *weights[next_id])

        partial_vals = new_partial

    result = np.empty(nodes.shape, order='F')
    for index in six.moves.xrange(degree + 1):
        key = (0,) * (degree - index) + (1,) * index
        result[index, :] = partial_vals[key]

    interval_delta = curve_end - curve_start
    true_start = curve_start + start * interval_delta
    true_end = curve_start + end * interval_delta
    return result, true_start, true_end
    # pylint: enable=too-many-locals


def _evaluate_hodograph(s, nodes):
    r"""Evaluate the Hodograph curve at a point :math:`s`.

    The Hodograph (first derivative) of a B |eacute| zier curve
    degree :math:`d = n - 1` and is given by

    .. math::

       B'(s) = n \sum_{j = 0}^{d} \binom{d}{j} s^j
       (1 - s)^{d - j} \cdot \Delta v_j

    where each forward difference is given by
    :math:`\Delta v_j = v_{j + 1} - v_j`.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes of a curve.
        s (float): A parameter along the curve at which the Hodograph
            is to be evaluated.

    Returns:
        numpy.ndarray: The point on the Hodograph curve (as a two
        dimensional NumPy array with a single row).
    """
    first_deriv = nodes[1:, :] - nodes[:-1, :]
    # NOTE: Taking the derivative drops the degree by 1.
    num_nodes, _ = nodes.shape
    degree = num_nodes - 1
    return degree * evaluate_multi(
        first_deriv, np.asfortranarray([s]))


def get_curvature(nodes, degree, tangent_vec, s):
    r"""Compute the signed curvature of a curve at :math:`s`.

    Computed via

    .. math::

       \frac{B'(s) \times B''(s)}{\left\lVert B'(s) \right\rVert_2^3}

    .. image:: images/get_curvature.png
       :align: center

    .. testsetup:: get-curvature

       import numpy as np
       import bezier
       from bezier._curve_helpers import evaluate_hodograph
       from bezier._curve_helpers import get_curvature

    .. doctest:: get-curvature
       :options: +NORMALIZE_WHITESPACE

       >>> nodes = np.asfortranarray([
       ...     [1.0 ,  0.0],
       ...     [0.75,  2.0],
       ...     [0.5 , -2.0],
       ...     [0.25,  2.0],
       ...     [0.0 ,  0.0],
       ... ])
       >>> s = 0.5
       >>> tangent_vec = evaluate_hodograph(s, nodes)
       >>> tangent_vec
       array([[-1., 0.]])
       >>> curvature = get_curvature(nodes, 4, tangent_vec, s)
       >>> curvature
       -12.0

    .. testcleanup:: get-curvature

       import make_images
       make_images.get_curvature(nodes, s, tangent_vec, curvature)

    Args:
        nodes (numpy.ndarray): The nodes of a curve.
        degree (int): The degree of the curve.
        tangent_vec (numpy.ndarray): The already computed value of
            :math:`B'(s)`
        s (float): The parameter value along the curve.

    Returns:
        float: The signed curvature.
    """
    if degree == 1:
        return 0.0

    # NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    #       here. It may be worthwhile to implement ``Curve.hodograph()``
    #       and ``Curve.concavity()`` to avoid re-computing the
    #       first and second node differences.
    first_deriv = nodes[1:, :] - nodes[:-1, :]
    second_deriv = first_deriv[1:, :] - first_deriv[:-1, :]
    concavity = degree * (degree - 1) * evaluate_multi(
        second_deriv, np.asfortranarray([s]))

    curvature = _helpers.cross_product(tangent_vec, concavity)
    # NOTE: We convert to 1D to make sure NumPy uses vector norm.
    curvature /= np.linalg.norm(tangent_vec[0, :], ord=2)**3
    return curvature


def newton_refine(curve, point, s):
    r"""Refine a solution to :math:`B(s) = p` using Newton's method.

    Computes updates via

    .. math::

       \mathbf{0} \approx
           \left(B\left(s_{\ast}\right) - p\right) +
           B'\left(s_{\ast}\right) \Delta s

    For example, consider the curve

    .. math::

       B(s) =
           \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)^2 +
           \left[\begin{array}{c} 1 \\ 2 \end{array}\right]
               2 (1 - s) s +
           \left[\begin{array}{c} 3 \\ 1 \end{array}\right] s^2

    and the point :math:`B\left(\frac{1}{4}\right) =
    \frac{1}{16} \left[\begin{array}{c} 9 \\ 13 \end{array}\right]`.

    Starting from the **wrong** point :math:`s = \frac{3}{4}`, we have

    .. math::

       \begin{align*}
       p - B\left(\frac{1}{2}\right) &= -\frac{1}{2}
           \left[\begin{array}{c} 3 \\ 1 \end{array}\right] \\
       B'\left(\frac{1}{2}\right) &= \frac{1}{2}
           \left[\begin{array}{c} 7 \\ -1 \end{array}\right] \\
       \Longrightarrow \frac{1}{4} \left[\begin{array}{c c}
           7 & -1 \end{array}\right] \left[\begin{array}{c}
           7 \\ -1 \end{array}\right] \Delta s &= -\frac{1}{4}
           \left[\begin{array}{c c} 7 & -1 \end{array}\right]
           \left[\begin{array}{c} 3 \\ 1 \end{array}\right] \\
       \Longrightarrow \Delta s &= -\frac{2}{5}
       \end{align*}

    .. image:: images/newton_refine_curve.png
       :align: center

    .. testsetup:: newton-refine-curve, newton-refine-curve-cusp

       import numpy as np
       import bezier
       from bezier._curve_helpers import newton_refine

    .. doctest:: newton-refine-curve
       :options: +NORMALIZE_WHITESPACE

       >>> nodes = np.asfortranarray([
       ...     [0.0, 0.0],
       ...     [1.0, 2.0],
       ...     [3.0, 1.0],
       ... ])
       >>> curve = bezier.Curve(nodes, degree=2)
       >>> point = curve.evaluate(0.25)
       >>> point
       array([[ 0.5625, 0.8125]])
       >>> s = 0.75
       >>> new_s = newton_refine(curve, point, s)
       >>> 5 * (new_s - s)
       -2.0

    .. testcleanup:: newton-refine-curve

       import make_images
       make_images.newton_refine_curve(curve, point, s, new_s)

    On curves that are not "valid" (i.e. :math:`B(s)` is not
    injective with non-zero gradient), Newton's method may
    break down and converge linearly:

    .. image:: images/newton_refine_curve_cusp.png
       :align: center

    .. doctest:: newton-refine-curve-cusp
       :options: +NORMALIZE_WHITESPACE

       >>> nodes = np.asfortranarray([
       ...     [ 6.0, -3.0],
       ...     [-2.0,  3.0],
       ...     [-2.0, -3.0],
       ...     [ 6.0,  3.0],
       ... ])
       >>> curve = bezier.Curve(nodes, degree=3)
       >>> expected = 0.5
       >>> point = curve.evaluate(expected)
       >>> point
       array([[ 0., 0.]])
       >>> s_vals = [0.625, None, None, None, None, None]
       >>> np.log2(abs(expected - s_vals[0]))
       -3.0
       >>> s_vals[1] = newton_refine(curve, point, s_vals[0])
       >>> np.log2(abs(expected - s_vals[1]))
       -3.983...
       >>> s_vals[2] = newton_refine(curve, point, s_vals[1])
       >>> np.log2(abs(expected - s_vals[2]))
       -4.979...
       >>> s_vals[3] = newton_refine(curve, point, s_vals[2])
       >>> np.log2(abs(expected - s_vals[3]))
       -5.978...
       >>> s_vals[4] = newton_refine(curve, point, s_vals[3])
       >>> np.log2(abs(expected - s_vals[4]))
       -6.978...
       >>> s_vals[5] = newton_refine(curve, point, s_vals[4])
       >>> np.log2(abs(expected - s_vals[5]))
       -7.978...

    .. testcleanup:: newton-refine-curve-cusp

       import make_images
       make_images.newton_refine_curve_cusp(curve, s_vals)

    Due to round-off, the Newton process terminates with an error that is not
    close to machine precision :math:`\varepsilon` when :math:`\Delta s = 0`.

    .. testsetup:: newton-refine-curve-cusp-continued

       import numpy as np
       import bezier
       from bezier._curve_helpers import newton_refine

       nodes = np.asfortranarray([
           [6.0, -3.0],
           [-2.0, 3.0],
           [-2.0, -3.0],
           [6.0, 3.0],
       ])
       curve = bezier.Curve(nodes, degree=3)
       point = curve.evaluate(0.5)

    .. doctest:: newton-refine-curve-cusp-continued

       >>> s_vals = [0.625]
       >>> new_s = newton_refine(curve, point, s_vals[-1])
       >>> while new_s not in s_vals:
       ...     s_vals.append(new_s)
       ...     new_s = newton_refine(curve, point, s_vals[-1])
       ...
       >>> terminal_s = s_vals[-1]
       >>> terminal_s == newton_refine(curve, point, terminal_s)
       True
       >>> 2.0**(-31) <= abs(terminal_s - 0.5) <= 2.0**(-29)
       True

    Due to round-off near the cusp, the final error resembles
    :math:`\sqrt{\varepsilon}` rather than machine precision
    as expected.

    Args:
        curve (.Curve): The curve to refine a point on.
        point (numpy.ndarray): A point on ``curve``.
        s (float): An "almost" solution to :math:`B(s) = p`.

    Returns:
        float: The updated value :math:`s + \Delta s`.
    """
    # pylint: disable=protected-access
    nodes = curve._nodes
    # pylint: enable=protected-access
    pt_delta = point - curve.evaluate(s)
    derivative = evaluate_hodograph(s, nodes)
    # Each array is 1x2 (i.e. a row vector), we want the vector dot product.
    delta_s = (
        np.vdot(pt_delta[0, :], derivative[0, :]) /
        np.vdot(derivative[0, :], derivative[0, :]))
    return s + delta_s


def locate_point(curve, point):
    r"""Locate a point on a curve.

    Does so by recursively subdividing the curve and rejecting
    sub-curves with bounding boxes that don't contain the point.
    After the sub-curves are sufficiently small, uses Newton's
    method to zoom in on the parameter value.

    .. note::

       This assumes, but does not check, that ``point`` is ``1xD``,
       where ``D`` is the dimension that ``curve`` is in.

    Args:
        curve (.Curve): A B |eacute| zier curve.
        point (numpy.ndarray): The point to locate.

    Returns:
        Optional[float]: The parameter value (:math:`s`) corresponding
        to ``point`` or :data:`None` if the point is not on the ``curve``.

    Raises:
        ValueError: If the standard deviation of the remaining start / end
            parameters among the subdivided intervals exceeds a given
            threshold (e.g. :math:`2^{-20}`).
    """
    candidates = [curve]
    for _ in six.moves.xrange(_MAX_LOCATE_SUBDIVISIONS + 1):
        next_candidates = []
        for candidate in candidates:
            nodes = candidate._nodes  # pylint: disable=protected-access
            if _helpers.contains_nd(nodes, point):
                next_candidates.extend(candidate.subdivide())

        candidates = next_candidates

    if not candidates:
        return None

    # pylint: disable=protected-access
    params = [
        (candidate._start, candidate._end) for candidate in candidates]
    # pylint: enable=protected-access

    if np.std(params) > _LOCATE_STD_CAP:
        raise ValueError(
            'Parameters not close enough to one another', params)

    s_approx = np.mean(params)
    return newton_refine(curve, point, s_approx)


def reduce_pseudo_inverse(nodes, degree):
    """Performs degree-reduction for a B |eacute| zier curve.

    Does so by using the pseudo-inverse of the degree elevation
    operator (which is overdetermined).

    Args:
        nodes (numpy.ndarray): The nodes in the curve.
        degree (int): The degree of the curve.

    Returns:
        numpy.ndarray: The reduced nodes.

    Raises:
        NotImplementedError: If the degree is not 1, 2, 3 or 4.
    """
    if degree == 1:
        reduction = _REDUCTION0
        denom = _REDUCTION_DENOM0
    elif degree == 2:
        reduction = _REDUCTION1
        denom = _REDUCTION_DENOM1
    elif degree == 3:
        reduction = _REDUCTION2
        denom = _REDUCTION_DENOM2
    elif degree == 4:
        reduction = _REDUCTION3
        denom = _REDUCTION_DENOM3
    else:
        raise NotImplementedError(degree)

    result = _helpers.matrix_product(reduction, nodes)
    result /= denom
    return result


def _projection_error(nodes, projected):
    """Compute the error between ``nodes`` and the projected nodes.

    Helper for :func:`_maybe_reduce`.

    For now, just compute the relative error in the Frobenius norm. But,
    we may wish to consider the error per row / point instead.

    Args:
        nodes (numpy.ndarray): Nodes in a curve.
        projected (numpy.ndarray): The ``nodes`` projected into the
            space of degree-elevated nodes.

    Returns:
        float: The relative error.
    """
    relative_err = np.linalg.norm(nodes - projected, ord='fro')
    if relative_err != 0.0:
        relative_err /= np.linalg.norm(nodes, ord='fro')

    return relative_err


def _maybe_reduce(nodes):
    r"""Reduce nodes in a curve if they are degree-elevated.

    We check if the nodes are degree-elevated by projecting onto the
    space of degree-elevated curves of the same degree, then comparing
    to the projection. We form the projection by taking the corresponding
    elevation matrix :math:`E` (from one degree lower) and forming
    :math:`E \left(E^T E\right)^{-1} E^T`.

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        Tuple[bool, numpy.ndarray]: Pair of values. The first indicates
        if the ``nodes`` were reduced. The second is the resulting nodes,
        either the reduced ones or the original passed in.

    Raises:
        NotImplementedError: If the curve is degree 5 or higher.
    """
    num_nodes, _ = nodes.shape
    if num_nodes < 2:
        return False, nodes
    elif num_nodes == 2:
        projection = _PROJECTION0
        denom = _PROJ_DENOM0
    elif num_nodes == 3:
        projection = _PROJECTION1
        denom = _PROJ_DENOM1
    elif num_nodes == 4:
        projection = _PROJECTION2
        denom = _PROJ_DENOM2
    elif num_nodes == 5:
        projection = _PROJECTION3
        denom = _PROJ_DENOM3
    else:
        raise NotImplementedError(num_nodes)

    projected = _helpers.matrix_product(projection, nodes) / denom
    relative_err = _projection_error(nodes, projected)
    if relative_err < _REDUCE_THRESHOLD:
        return True, reduce_pseudo_inverse(nodes, num_nodes - 1)
    else:
        return False, nodes


def full_reduce(nodes):
    """Apply degree reduction to ``nodes`` until it can no longer be reduced.

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        numpy.ndarray: The fully degree-reduced nodes.
    """
    was_reduced, nodes = _maybe_reduce(nodes)
    while was_reduced:
        was_reduced, nodes = _maybe_reduce(nodes)
    return nodes


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    evaluate_multi_barycentric = _evaluate_multi_barycentric
    evaluate_multi = _evaluate_multi
    specialize_curve = _specialize_curve
    evaluate_hodograph = _evaluate_hodograph
else:
    evaluate_multi_barycentric = _speedup.speedup.evaluate_curve_barycentric
    evaluate_multi = _speedup.speedup.evaluate_multi
    specialize_curve = _speedup.speedup.specialize_curve
    evaluate_hodograph = _speedup.speedup.evaluate_hodograph
# pylint: enable=invalid-name
