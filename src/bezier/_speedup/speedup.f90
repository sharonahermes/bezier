! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     https://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module speedup

  use iso_c_binding, only: c_double, c_int, c_bool
  implicit none
  private
  public &
       de_casteljau_one_round, evaluate_curve_barycentric, evaluate_multi, &
       linearization_error, evaluate_barycentric, evaluate_barycentric_multi, &
       evaluate_cartesian_multi, cross_product, segment_intersection, bbox, &
       specialize_curve_generic, specialize_curve_quadratic, &
       specialize_curve, jacobian_both, evaluate_hodograph, &
       newton_refine_intersect, jacobian_det, bbox_intersect, &
       wiggle_interval, parallel_different, from_linearized

  integer, parameter :: dp=kind(0.d0)

contains

  subroutine de_casteljau_one_round( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, new_nodes) &
       bind(c, name='de_casteljau_one_round')

    ! NOTE: This is de Casteljau on a Bezier surface / triangle.

    integer(c_int), intent(in), value :: num_nodes
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in), value :: degree
    real(c_double), intent(in), value :: lambda1
    real(c_double), intent(in), value :: lambda2
    real(c_double), intent(in), value :: lambda3
    real(c_double), intent(out) :: new_nodes(num_nodes - degree - 1, dimension_)
    ! Variables outside of signature.
    integer :: index_
    integer :: parent_i1
    integer :: parent_i2
    integer :: parent_i3
    integer :: k, j

    index_ = 1
    parent_i1 = 1
    parent_i2 = 2
    parent_i3 = degree + 2
    ! NOTE: Throughout for index_ <--> (i, j, k) in the (degree - 1)
    !       triangle, we have parent_i1 = index_ + k <--> (i + 1, j, k),
    !       parent_i2 = index_ + k + 1 <--> (i, j + 1, k) and
    !       parent_i3 = index_ + degree + 1 <--> (i, j, k + 1).
    do k = 0, degree - 1
       do j = 0, degree - k - 1
          ! NOTE: i = (degree - 1) - j - k
          new_nodes(index_, :) = ( &
               lambda1 * nodes(parent_i1, :) + &
               lambda2 * nodes(parent_i2, :) + &
               lambda3 * nodes(parent_i3, :))
          ! Update all the indices.
          parent_i1 = parent_i1 + 1
          parent_i2 = parent_i2 + 1
          parent_i3 = parent_i3 + 1
          index_ = index_ + 1
       end do
       ! Update the indices that depend on k.
       parent_i1 = parent_i1 + 1
       parent_i2 = parent_i2 + 1
    end do

  end subroutine de_casteljau_one_round

  subroutine evaluate_curve_barycentric( &
       degree, dimension_, nodes, num_vals, lambda1, lambda2, evaluated) &
       bind(c, name='evaluate_curve_barycentric')

    ! NOTE: This is evaluate_multi_barycentric for a Bezier curve.

    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(degree + 1, dimension_)
    integer(c_int), intent(in), value :: num_vals
    real(c_double), intent(in) :: lambda1(num_vals)
    real(c_double), intent(in) :: lambda2(num_vals)
    real(c_double), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer :: i, j
    real(dp) :: lambda2_pow(num_vals)
    integer :: binom_val

    lambda2_pow = 1.0_dp
    binom_val = 1

    forall (i = 1:num_vals)
       evaluated(i, :) = lambda1(i) * nodes(1, :)
    end forall

    do i = 2, degree
       lambda2_pow = lambda2_pow * lambda2
       binom_val = (binom_val * (degree - i + 2)) / (i - 1)
       forall (j = 1:num_vals)
          evaluated(j, :) = ( &
               evaluated(j, :) + &
               binom_val * lambda2_pow(j) * nodes(i, :)) * lambda1(j)
       end forall
    end do

    forall (i = 1:num_vals)
       evaluated(i, :) = ( &
            evaluated(i, :) + &
            lambda2_pow(i) * lambda2(i) * nodes(degree + 1, :))
    end forall

  end subroutine evaluate_curve_barycentric

  subroutine evaluate_multi( &
       degree, dimension_, nodes, num_vals, s_vals, evaluated) &
       bind(c, name='evaluate_multi')

    ! NOTE: This is evaluate_multi for a Bezier curve.

    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(degree + 1, dimension_)
    integer(c_int), intent(in), value :: num_vals
    real(c_double), intent(in) :: s_vals(num_vals)
    real(c_double), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    real(dp) :: one_less(num_vals)

    one_less = 1.0_dp - s_vals
    call evaluate_curve_barycentric( &
         degree, dimension_, nodes, num_vals, one_less, s_vals, evaluated)
  end subroutine evaluate_multi

  subroutine linearization_error(degree, dimension_, nodes, error) &
       bind(c, name='linearization_error')

    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(degree + 1, dimension_)
    real(c_double), intent(out) :: error
    ! Variables outside of signature.
    real(dp) :: second_deriv(degree - 1, dimension_)
    real(dp) :: worst_case(dimension_)

    if (degree == 1) then
       error = 0.0_dp
       return
    end if

    second_deriv = ( &
         nodes(:degree - 1, :) - &
         2.0_dp * nodes(2:degree, :) + &
         nodes(3:, :))
    worst_case = maxval(abs(second_deriv), 1)
    error = 0.125_dp * degree * (degree - 1) * norm2(worst_case)
  end subroutine linearization_error

  subroutine evaluate_barycentric( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, point) &
       bind(c, name='evaluate_barycentric')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    integer(c_int), intent(in), value :: num_nodes
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in), value :: degree
    real(c_double), intent(in), value :: lambda1
    real(c_double), intent(in), value :: lambda2
    real(c_double), intent(in), value :: lambda3
    real(c_double), intent(out) :: point(1, dimension_)
    ! Variables outside of signature.
    real(dp) :: param_vals(1, 3)

    param_vals(1, 1) = lambda1
    param_vals(1, 2) = lambda2
    param_vals(1, 3) = lambda3
    call evaluate_barycentric_multi( &
         num_nodes, dimension_, nodes, degree, 1, param_vals, point)

  end subroutine evaluate_barycentric

  subroutine evaluate_barycentric_multi( &
       num_nodes, dimension_, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='evaluate_barycentric_multi')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    integer(c_int), intent(in), value :: num_nodes
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 3)
    real(c_double), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer :: k, binom_val, index_, new_index
    real(dp) :: row_result(num_vals, dimension_)

    index_ = num_nodes
    forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
       evaluated(new_index, :) = nodes(index_, :)
    end forall

    if (degree == 0) then
       return
    end if

    binom_val = 1
    do k = degree - 1, 0, -1
        ! We want to go from (d C (k + 1)) to (d C k).
        binom_val = (binom_val * (k + 1)) / (degree - k)
        index_ = index_ - 1  ! Step to last element in row.
        !     k = d - 1, d - 2, ...
        ! d - k =     1,     2, ...
        ! We know row k has (d - k + 1) elements.
        new_index = index_ - degree + k  ! First element in row.

        ! lambda1 = param_vals(:, 1)
        ! lambda2 = param_vals(:, 2)
        call evaluate_curve_barycentric( &
             degree - k, dimension_, nodes(new_index:index_, :), &
             num_vals, param_vals(:, 1), param_vals(:, 2), row_result)

        ! Update index for next iteration.
        index_ = new_index

        ! lambda3 = param_vals(:, 3)
        forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
           evaluated(new_index, :) = ( &
                param_vals(new_index, 3) * evaluated(new_index, :) + &
                binom_val * row_result(new_index, :))
        end forall
    end do

  end subroutine evaluate_barycentric_multi

  subroutine evaluate_cartesian_multi( &
       num_nodes, dimension_, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='evaluate_cartesian_multi')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This mostly copies evaluate_barycentric_multi but does not just
    !       call it directly. This is to avoid copying param_vals.

    integer(c_int), intent(in), value :: num_nodes
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 2)
    real(c_double), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer :: k, binom_val, index_, new_index
    real(dp) :: row_result(num_vals, dimension_)
    real(dp) :: lambda1_vals(num_vals)

    index_ = num_nodes
    forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
       evaluated(new_index, :) = nodes(index_, :)
    end forall

    if (degree == 0) then
       return
    end if

    lambda1_vals = 1.0_dp - param_vals(:, 1) - param_vals(:, 2)

    binom_val = 1
    do k = degree - 1, 0, -1
        ! We want to go from (d C (k + 1)) to (d C k).
        binom_val = (binom_val * (k + 1)) / (degree - k)
        index_ = index_ - 1  ! Step to last element in row.
        !     k = d - 1, d - 2, ...
        ! d - k =     1,     2, ...
        ! We know row k has (d - k + 1) elements.
        new_index = index_ - degree + k  ! First element in row.

        ! lambda1 = param_vals(:, 1)
        ! lambda2 = param_vals(:, 1)
        call evaluate_curve_barycentric( &
             degree - k, dimension_, nodes(new_index:index_, :), &
             num_vals, lambda1_vals, param_vals(:, 1), row_result)

        ! Update index for next iteration.
        index_ = new_index

        ! lambda3 = param_vals(:, 2)
        forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
           evaluated(new_index, :) = ( &
                param_vals(new_index, 2) * evaluated(new_index, :) + &
                binom_val * row_result(new_index, :))
        end forall
    end do

  end subroutine evaluate_cartesian_multi

  subroutine cross_product(vec0, vec1, result_) bind(c, name='cross_product')

    real(c_double), intent(in) :: vec0(1, 2)
    real(c_double), intent(in) :: vec1(1, 2)
    real(c_double), intent(out) :: result_

    result_ = vec0(1, 1) * vec1(1, 2) - vec0(1, 2) * vec1(1, 1)

  end subroutine cross_product

  subroutine segment_intersection(start0, end0, start1, end1, s, t, success) &
       bind(c, name='segment_intersection')

    real(c_double), intent(in) :: start0(1, 2)
    real(c_double), intent(in) :: end0(1, 2)
    real(c_double), intent(in) :: start1(1, 2)
    real(c_double), intent(in) :: end1(1, 2)
    real(c_double), intent(out) :: s, t
    logical(c_bool), intent(out) :: success
    ! Variables outside of signature.
    real(dp) :: delta0(1, 2)
    real(dp) :: delta1(1, 2)
    real(dp) :: start_delta(1, 2)
    real(dp) :: cross_d0_d1
    real(dp) :: other_cross

    delta0 = end0 - start0
    delta1 = end1 - start1
    call cross_product(delta0, delta1, cross_d0_d1)

    if (cross_d0_d1 == 0.0_dp) then
       success = .FALSE.
    else
       start_delta = start1 - start0
       call cross_product(start_delta, delta1, other_cross)
       s = other_cross / cross_d0_d1
       call cross_product(start_delta, delta0, other_cross)
       t = other_cross / cross_d0_d1
       success = .TRUE.
    end if

  end subroutine segment_intersection

  subroutine bbox(num_nodes, nodes, left, right, bottom, top) &
       bind(c, name='bbox')

    integer(c_int), intent(in), value :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    real(c_double), intent(out) :: left, right, bottom, top
    ! Variables outside of signature.
    real(dp) :: workspace(2)

    workspace = minval(nodes, 1)
    left = workspace(1)
    bottom = workspace(2)
    workspace = maxval(nodes, 1)
    right = workspace(1)
    top = workspace(2)

  end subroutine bbox

  subroutine specialize_curve_generic( &
       nodes, degree, dimension_, start, end_, new_nodes)

    ! NOTE: This is a helper for ``specialize_curve`` that works on any degree.

    real(dp), intent(in) :: nodes(degree + 1, dimension_)
    integer :: dimension_
    integer, intent(in) :: degree
    real(dp), intent(in) :: start, end_
    real(dp), intent(out) :: new_nodes(degree + 1, dimension_)
    ! Variables outside of signature.
    real(dp) :: workspace(degree, dimension_, degree + 1)
    integer :: index_, curr_size, j
    real(dp) :: minus_start, minus_end

    minus_start = 1.0_dp - start
    minus_end = 1.0_dp - end_
    workspace(:, :, 1) = minus_start * nodes(:degree, :) + start * nodes(2:, :)
    workspace(:, :, 2) = minus_end * nodes(:degree, :) + end_ * nodes(2:, :)

    curr_size = degree
    do index_ = 3, degree + 1
       curr_size = curr_size - 1
       ! First add a new "column" (or whatever the 3rd dimension is called)
       ! at the end using ``end_``.
       workspace(:curr_size, :, index_) = ( &
            minus_end * workspace(:curr_size, :, index_ - 1) + &
            end_ * workspace(2:curr_size + 1, :, index_ - 1))
       ! Update all the values in place by using de Casteljau with the
       ! ``start`` parameter.
       forall (j = 1:index_ - 1)
          workspace(:curr_size, :, j) = ( &
               minus_start * workspace(:curr_size, :, j) + &
               start * workspace(2:curr_size + 1, :, j))
       end forall
    end do

    ! Move the final "column" (or whatever the 3rd dimension is called)
    ! of the workspace into ``new_nodes``.
    forall (index_ = 1:degree + 1)
       new_nodes(index_, :) = workspace(1, :, index_)
    end forall

  end subroutine specialize_curve_generic

  subroutine specialize_curve_quadratic( &
       nodes, dimension_, start, end_, new_nodes)

    real(dp), intent(in) :: nodes(3, dimension_)
    integer :: dimension_
    real(dp), intent(in) :: start, end_
    real(dp), intent(out) :: new_nodes(3, dimension_)
    ! Variables outside of signature.
    real(dp) :: minus_start, minus_end, prod_both

    minus_start = 1.0_dp - start
    minus_end = 1.0_dp - end_
    prod_both = start * end_

    new_nodes(1, :) = ( &
         minus_start * minus_start * nodes(1, :) + &
         2.0_dp * start * minus_start * nodes(2, :) + &
         start * start * nodes(3, :))
    new_nodes(2, :) = ( &
         minus_start * minus_end * nodes(1, :) + &
         (end_ + start - 2.0_dp * prod_both) * nodes(2, :) + &
         prod_both * nodes(3, :))
    new_nodes(3, :) = ( &
         minus_end * minus_end * nodes(1, :) + &
         2.0_dp * end_ * minus_end * nodes(2, :) + &
         end_ * end_ * nodes(3, :))

  end subroutine specialize_curve_quadratic

  subroutine specialize_curve( &
       degree, dimension_, nodes, start, end_, curve_start, curve_end, &
       new_nodes, true_start, true_end) bind(c, name='specialize_curve')

    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(degree + 1, dimension_)
    real(c_double), intent(in), value :: start, end_, curve_start, curve_end
    real(c_double), intent(out) :: new_nodes(degree + 1, dimension_)
    real(c_double), intent(out) :: true_start, true_end
    ! Variables outside of signature.
    real(dp) :: interval_delta

    if (degree == 1) then
       new_nodes(1, :) = (1.0_dp - start) * nodes(1, :) + start * nodes(2, :)
       new_nodes(2, :) = (1.0_dp - end_) * nodes(1, :) + end_ * nodes(2, :)
    else if (degree == 2) then
       call specialize_curve_quadratic( &
            nodes, dimension_, start, end_, new_nodes)
    else
       call specialize_curve_generic( &
            nodes, degree, dimension_, start, end_, new_nodes)
    end if

    ! Now, compute the new interval.
    interval_delta = curve_end - curve_start
    true_start = curve_start + start * interval_delta
    true_end = curve_start + end_ * interval_delta

  end subroutine specialize_curve

  subroutine jacobian_both( &
       num_nodes, dimension_, nodes, degree, new_nodes) &
       bind(c, name='jacobian_both')

    integer(c_int), intent(in), value :: num_nodes
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in), value :: degree
    real(c_double), intent(out) :: &
         new_nodes(num_nodes - degree - 1, 2 * dimension_)
    ! Variables outside of signature.
    integer :: index_, i, j, k, num_vals

    index_ = 1
    i = 1
    j = degree + 2
    do num_vals = degree, 1, -1
       do k = 0, num_vals - 1
          ! jacobian_s
          new_nodes(index_, :dimension_) = nodes(i + 1, :) - nodes(i, :)
          ! jacobian_t
          new_nodes(index_, dimension_ + 1:) = nodes(j, :) - nodes(i, :)
          ! Update the indices
          index_ = index_ + 1
          i = i + 1
          j = j + 1
       end do
       ! In between each row, the index_ gains an extra value.
       i = i + 1
    end do

    new_nodes = degree * new_nodes

  end subroutine jacobian_both

  subroutine evaluate_hodograph(s, degree, dimension_, nodes, hodograph) &
       bind(c, name='evaluate_hodograph')

    real(c_double), intent(in), value :: s
    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: dimension_
    real(c_double), intent(in) :: nodes(degree + 1, dimension_)
    real(c_double), intent(out) :: hodograph(1, dimension_)
    ! Variables outside of signature.
    real(dp) :: first_deriv(degree, dimension_)
    real(dp) :: param(1)

    first_deriv = nodes(2:, :) - nodes(:degree, :)
    param = s
    call evaluate_multi( &
         degree - 1, dimension_, first_deriv, 1, param, hodograph)
    hodograph = degree * hodograph

  end subroutine evaluate_hodograph

  subroutine newton_refine_intersect( &
       s, degree1, nodes1, t, degree2, nodes2, new_s, new_t) &
       bind(c, name='newton_refine_intersect')

    real(c_double), intent(in), value :: s
    integer(c_int), intent(in), value :: degree1
    real(c_double), intent(in) :: nodes1(degree1 + 1, 2)
    real(c_double), intent(in), value :: t
    integer(c_int), intent(in), value :: degree2
    real(c_double), intent(in) :: nodes2(degree2 + 1, 2)
    real(c_double), intent(out) :: new_s, new_t
    ! Variables outside of signature.
    real(dp) :: param(1)
    real(dp) :: func_val(1, 2)
    real(dp) :: workspace(1, 2)
    real(dp) :: jac_mat(2, 2)
    real(dp) :: determinant, delta_s, delta_t

    param = t
    call evaluate_multi( &
         degree2, 2, nodes2, 1, param, func_val)
    param = s
    call evaluate_multi( &
         degree1, 2, nodes1, 1, param, workspace)
    func_val = func_val - workspace

    if (all(func_val == 0.0_dp)) then
       new_s = s
       new_t = t
       return
    end if

    call evaluate_hodograph(s, degree1, 2, nodes1, jac_mat(1:1, :))
    ! NOTE: We actually want the negative, since we want -B2'(t), but
    !       since we manually solve the system, it's just algebra
    !       to figure out how to use the negative values.
    call evaluate_hodograph(t, degree2, 2, nodes2, jac_mat(2:2, :))

    determinant = ( &
         jac_mat(1, 1) * jac_mat(2, 2) - jac_mat(1, 2) * jac_mat(2, 1))

    ! NOTE: We manually invert the 2x2 system ([ds, dt] J)^T = f^T.
    delta_s = ( &
         (jac_mat(2, 2) * func_val(1, 1) - &
         jac_mat(2, 1) * func_val(1, 2)) / determinant)
    new_s = s + delta_s

    delta_t = ( &
         (jac_mat(1, 2) * func_val(1, 1) - &
         jac_mat(1, 1) * func_val(1, 2)) / determinant)
    new_t = t + delta_t

  end subroutine newton_refine_intersect

  subroutine jacobian_det( &
       num_nodes, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='jacobian_det')

    integer(c_int), intent(in), value :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    integer(c_int), intent(in), value :: degree
    integer(c_int), intent(in), value :: num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 2)
    real(c_double), intent(out) :: evaluated(num_vals)
    ! Variables outside of signature.
    real(dp) :: jac_nodes(num_nodes - degree - 1, 4)
    real(dp) :: Bs_Bt_vals(num_vals, 4)
    real(dp) :: determinant

    call jacobian_both( &
         num_nodes, 2, nodes, degree, jac_nodes)
    if (degree == 1) then
       determinant = ( &
            jac_nodes(1, 1) * jac_nodes(1, 4) - &
            jac_nodes(1, 2) * jac_nodes(1, 3))
       evaluated = determinant
    else
       call evaluate_cartesian_multi( &
            num_nodes - degree - 1, 4, jac_nodes, degree - 1, &
            num_vals, param_vals, Bs_Bt_vals)
       evaluated = ( &
            Bs_Bt_vals(:, 1) * Bs_Bt_vals(:, 4) - &
            Bs_Bt_vals(:, 2) * Bs_Bt_vals(:, 3))
    end if

  end subroutine jacobian_det

  subroutine bbox_intersect(num_nodes1, nodes1, num_nodes2, nodes2, enum_) &
       bind(c, name='bbox_intersect')

    integer(c_int), intent(in), value :: num_nodes1, num_nodes2
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    integer(c_int), intent(out) :: enum_
    ! Variables outside of signature.
    real(dp) :: left1, right1, bottom1, top1
    real(dp) :: left2, right2, bottom2, top2

    call bbox(num_nodes1, nodes1, left1, right1, bottom1, top1)
    call bbox(num_nodes2, nodes2, left2, right2, bottom2, top2)

    if ( &
         right2 < left1 .OR. right1 < left2 .OR. &
         top2 < bottom1 .OR. top1 < bottom2) then
        enum_ = 2
     else if ( &
          right2 == left1 .OR. right1 == left2 .OR. &
          top2 == bottom1 .OR. top1 == bottom2) then
        enum_ = 1
    else
        enum_ = 0
    end if

  end subroutine bbox_intersect

  subroutine wiggle_interval(value_, result_, success) &
       bind(c, name='wiggle_interval')

    real(c_double), intent(in), value :: value_
    real(c_double), intent(out) :: result_
    logical(c_bool), intent(out) :: success

    success = .TRUE.
    if (-0.5_dp**45 < value_ .AND. value_ < 0.5_dp**45) then
       result_ = 0.0_dp
    else if (0.5_dp**45 <= value_ .AND. value_ <= 1.0_dp - 0.5_dp**45) then
       result_ = value_
    else if ( &
         1.0_dp - 0.5_dp**45 < value_ .AND. value_ < 1.0_dp + 0.5_dp**45) then
       result_ = 1.0_dp
    else
       success = .FALSE.
    end if

  end subroutine wiggle_interval

  subroutine parallel_different(start0, end0, start1, end1, result_) &
       bind(c, name='parallel_different')

    real(c_double), intent(in) :: start0(1, 2)
    real(c_double), intent(in) :: end0(1, 2)
    real(c_double), intent(in) :: start1(1, 2)
    real(c_double), intent(in) :: end1(1, 2)
    logical(c_bool), intent(out) :: result_
    ! Variables outside of signature.
    real(dp) :: delta0(1, 2)
    real(dp) :: val1, val2, val3

    delta0 = end0 - start0
    call cross_product(start0, delta0, val1)  ! line0_const
    call cross_product(start1, delta0, val2)  ! start1_against

    if (val1 /= val2) then
       result_ = .TRUE.
       return
    end if

    val1 = dot_product(delta0(1, :), delta0(1, :))  ! norm0_sq
    val2 = dot_product(start1(1, :) - start0(1, :), delta0(1, :))  ! start_numer
    !      0 <= start_numer / norm0_sq <= 1
    ! <==> 0 <= start_numer            <= norm0_sq
    if (0.0_dp <= val2 .AND. val2 <= val1) then
       result_ = .FALSE.
       return
    end if

    val2 = dot_product(end1(1, :) - start0(1, :), delta0(1, :))  ! end_numer
    !      0 <= end_numer / norm0_sq <= 1
    ! <==> 0 <= end_numer            <= norm0_sq
    if (0.0_dp <= val2 .AND. val2 <= val1) then
       result_ = .FALSE.
       return
    end if

    ! We know neither the start or end parameters are in [0, 1], but
    ! they may contain [0, 1] between them.
    val3 = min(val1, val2)  ! min_val
    val1 = max(val1, val2)  ! max_val

    ! So we make sure that 0 isn't between them.
    if (val3 <= 0.0_dp .AND. 0.0_dp <= val1) then
       result_ = .FALSE.
    else
       result_ = .TRUE.
    end if

  end subroutine parallel_different

  subroutine from_linearized( &
       error1, start1, end1, start_node1, end_node1, degree1, nodes1, &
       error2, start2, end2, start_node2, end_node2, degree2, nodes2, &
       refined_s, refined_t, does_intersect, py_exc) &
       bind(c, name='from_linearized')

    real(c_double), intent(in), value :: error1, start1, end1
    real(c_double), intent(in) :: start_node1(1, 2)
    real(c_double), intent(in) :: end_node1(1, 2)
    integer(c_int), intent(in), value :: degree1
    real(c_double), intent(in) :: nodes1(degree1 + 1, 2)
    real(c_double), intent(in), value :: error2, start2, end2
    real(c_double), intent(in) :: start_node2(1, 2)
    real(c_double), intent(in) :: end_node2(1, 2)
    integer(c_int), intent(in), value :: degree2
    real(c_double), intent(in) :: nodes2(degree2 + 1, 2)
    real(c_double), intent(out) :: refined_s, refined_t
    logical(c_bool), intent(out) :: does_intersect
    integer(c_int), intent(out) :: py_exc
    ! Variables outside of signature.
    real(dp) :: s, t
    integer :: enum_
    logical(1) :: success

    py_exc = 0
    call segment_intersection( &
         start_node1, end_node1, start_node2, end_node2, s, t, success)

    if (success) then
       ! Special case for lines, allow no leeway on almost intersections.
       if (error1 == 0.0_dp .AND. (s < 0.0_dp .OR. 1.0_dp < s)) then
          does_intersect = .FALSE.
          return
       end if

       if (error2 == 0.0_dp .AND. (t < 0.0_dp .OR. 1.0_dp < t)) then
          does_intersect = .FALSE.
          return
       end if

       if (s < -(0.5_dp**16) .OR. 1.0_dp + 0.5_dp**16 < s) then
          does_intersect = .FALSE.
          return
       end if

       if (t < -(0.5_dp**16) .OR. 1.0_dp + 0.5_dp**16 < t) then
          does_intersect = .FALSE.
          return
       end if
    else
       ! Handle special case where the curves are actually lines.
       if (error1 == 0.0_dp .AND. error2 == 0.0_dp) then
          call parallel_different( &
               start_node1, end_node1, start_node2, end_node2, success)
          if (success) then
             does_intersect = .FALSE.
             return
          end if
       else
          call bbox_intersect( &
               degree1 + 1, nodes1, degree2 + 1, nodes2, enum_)
          if (enum_ == 2) then  ! Disjoint
             does_intersect = .FALSE.
             return
          end if
       end if

       ! Expect the wrapper code to raise
       ! NotImplementedError('Line segments parallel.')
       py_exc = 1
       return
    end if

    does_intersect = .TRUE.
    ! Now, promote ``s`` and ``t`` onto the original curves.
    s = (1.0_dp - s) * start1 + s * end1  ! orig_s
    t = (1.0_dp - t) * start2 + t * end2  ! orig_t
    ! Perform one step of Newton iteration to refine the computed
    ! values of s and t.
    call newton_refine_intersect( &
         s, degree1, nodes1, t, degree2, nodes2, refined_s, refined_t)

    call wiggle_interval(refined_s, s, success)
    if (.NOT. success) then
       ! py_exc==2 indicates ``wiggle_interval`` failed.
       py_exc = 2
       return
    end if
    refined_s = s

    call wiggle_interval(refined_t, t, success)
    if (.NOT. success) then
       ! py_exc==2 indicates ``wiggle_interval`` failed.
       py_exc = 2
       return
    end if
    refined_t = t

  end subroutine from_linearized

end module speedup
