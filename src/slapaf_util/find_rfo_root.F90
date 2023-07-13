!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2013,2018, Ignacio Fdez. Galvan                        *
!***********************************************************************
!  Find_RFO_Root
!
!> @brief
!>   Find the \f$ \alpha \f$ value for the RS-RFO family of algorithms.
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Find the \f$ \alpha \f$ value that produces a given step length in the
!> RS-RFO family of algorithms. \cite Bes1998-TCA-100-265
!> This routine assumes that a larger \f$ \alpha \f$ produces a smaller step.
!>
!> Use a *regula falsi*-like method, instead of diagonalization as in eqs.
!> (16) and (20). This subroutine is designed to be used in an iterative loop.
!>
!> (\p x1, \p y1) and (\p x2, \p y2) are two \f$ \alpha \f$ and step length
!> values bracketing the desired root. If \p x2 = `0.0`, an appropriate value
!> for \p x2 will first be searched. (\p x3, \p y3) is another point between
!> \p x1 and \p x2. On output a new value for \p x3 will be suggested (and
!> the bracketing points will be updated), the surrounding loop must calculate
!> the new \p y3 value.
!>
!> @param[in,out] x1,y1     Smaller \f$ \alpha \f$ and larger step length
!> @param[in,out] x2,y2     Larger \f$ \alpha \f$ and smaller step length
!> @param[in,out] x3        Middle \f$ \alpha \f$, new value on output
!> @param[in]     y3        Step length at \p x3
!> @param[in]     Val       Desired step length value
!>
!> @see ::RS_RFO
!> @see ::RS_I_RFO
!> @see ::RS_P_RFO
!***********************************************************************

subroutine Find_RFO_Root(x1,y1,x2,y2,x3,y3,Val)

use Constants, only: Zero, One, Two, Four, Half, OneHalf
use Definitions, only: wp

implicit none
real(kind=wp), intent(inout) :: x1, y1, x2, y2, x3
real(kind=wp), intent(in) :: y3, Val
real(kind=wp) :: coefA, coefB, coefC, delta, denom, discr, new, nnew, xx1, xx2, yy1, yy2
real(kind=wp), parameter :: Thr = 1.0e-16_wp

! First, we must find a value of alpha that gives a small enough step

if (y2 > Val) then
  y2 = y3

  ! In the first iteration, simply try with alpha+1
  if (x2 == Zero) then
    new = x1+One
    x2 = new

    ! If the step is below the threshold, a bracketing pair has been found.
    ! Suggest an intermediate point (with secant method, or midpoint)
  else if (y3 < Val) then
    new = x1+(Val-y1)/(y2-y1)*(x2-x1)
    if ((new <= x1) .or. (new >= x2)) new = Half*(x1+x2)

    ! If it is not there yet, extrapolate using a straight line (times an
    ! arbitrary factor of 1.5, to overcome systematic undershooting), and
    ! update the other end point
  else
    if (y1-y2 > Thr) then
      ! Don't get too ambitious -- thus the min function.
      delta = min(x2,(Val-y2)*(x1-x2)/(y1-y2))
      new = x2+OneHalf*delta
    else
      ! If y2 is very close to y1, or actually lower, just double the step.
      new = x2+Two*(x2-x1)
    end if
    x1 = x2
    y1 = y2
    x2 = new
  end if

else

  ! Once a bracketing pair has been found, we arrive here with two end points
  ! (long and short), and a middle point, for which the step length has just
  ! been computed

  ! Use the midpoint (bisection) or secant (regula falsi) as safety net
  ! Do this in the subinterval where the root should be
  xx1 = x1
  yy1 = y1
  xx2 = x2
  yy2 = y2
  if (y3 < Val) then
    xx2 = x3
    yy2 = y3
  else
    xx1 = x3
    yy1 = y3
  end if
  ! If the middle point is shorter than the "short" step, something went
  ! wrong, start again from 1 (signal with -1),
  ! maybe now we will get a better eigenvalue to start with
  if ((y3 < y2) .and. (y3 > Val)) then
    x3 = -One
    return
  end if
  new = xx1+(Val-yy1)/(yy2-yy1)*(xx2-xx1)
  if ((new <= xx1) .or. (new >= xx2)) new = Half*(xx1+xx2)
  nnew = new

  ! Define a quadratic fitting for the 3 points
  denom = (x1-x2)*(x1-x3)*(x2-x3)
  if (abs(denom) > Thr) then
    coefA = (x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))/denom
    coefB = (x3**2*(y1-y2)+x2**2*(y3-y1)+x1**2*(y2-y3))/denom
    coefC = (x2*x3*(x2-x3)*y1+x3*x1*(x3-x1)*y2+x1*x2*(x1-x2)*y3)/denom
    discr = coefB**2-Four*coefA*(coefC-Val)
  else
    discr = Zero
  end if

  ! Find the point between x1 and x2 where the step would be exactly Val
  if (discr > Zero) then
    ! Choose root depending on sign of y1-y2
    if (y1-y2 > Zero) then
      nnew = (-coefB-sqrt(discr))/(Two*coefA)
    else if (y1-y2 < Zero) then
      nnew = (-coefB+sqrt(discr))/(Two*coefA)
    end if
  end if
  ! Last check: make sure the new point is in the correct interval
  if ((nnew > xx1) .and. (nnew < xx2)) new = nnew

  ! Update the end points
  x1 = xx1
  y1 = yy1
  x2 = xx2
  y2 = yy2
end if

! Update the suggested value
x3 = new

return

end subroutine Find_RFO_Root
