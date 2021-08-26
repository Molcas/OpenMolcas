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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine MnBrak2(ax,bx,cx,fa,fb,fc,f,q_a,q_b,dipole_a,dipole_b,r_a,r_b)

use Constants, only: Zero, One, Five, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: ax, bx
real(kind=wp), intent(out) :: cx, fa, fb, fc
! External function f and its arguments
real(kind=wp), external :: f
real(kind=wp), intent(in) :: q_a, q_b, dipole_a, dipole_b, r_a, r_b
real(kind=wp) :: vx, fv, coefA, coefB
logical(kind=iwp) :: Def
real(kind=wp), parameter :: Lim = 100.0_wp, Ratio = Half*(One+sqrt(Five)), Thr = 1.0e-20_wp
logical(kind=iwp), parameter :: Absolute = .true.

fa = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,ax,Absolute)
fb = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,bx,Absolute)
if (fa < fb) then
  cx = ax
  ax = bx
  bx = cx
  fc = fa
  fa = fb
  fb = fc
end if
! three points such that b is between a and c,
! and f(a) > f(b) > f(c), stop when f(c) > f(b)
cx = bx+Ratio*(bx-ax)
fc = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,cx,Absolute)
do while (fc <= fb)
  write(u6,*) ax,bx,cx
  Def = .true.
  ! try a parabolic fitting
  coefA = (cx*(fb-fa)+bx*(fa-fc)+ax*(fc-fb))
  coefB = (cx**2*(fa-fb)+bx**2*(fc-fa)+ax**2*(fb-fc))
  ! only worth it if the points are not linear and
  ! if the 2nd derivative is positive
  if ((abs(coefA) > Thr) .and. (coefA*(ax-cx) > Zero)) then
    vx = -Half*coefB/coefA
    ! v is between b and c
    if ((cx-vx)*(vx-bx) > Zero) then
      fv = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,vx,Absolute)
      ! minimum between b and c
      if (fv < fc) then
        ax = bx
        bx = vx
        fa = fb
        fb = fv
        return
        ! minimum between a and v
      else if (fv > fb) then
        cx = vx
        fc = fv
        return
      end if
      ! v is beyond c, but within limits
    else if ((bx+Lim*(cx-bx)-vx)*(vx-cx) > Zero) then
      fv = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,vx,Absolute)
      ! whatever happens, replace c,v -> b,c
      bx = cx
      cx = vx
      fb = fc
      fc = fv
      ! minimum between b and v
      if (fv > fc) then
        ax = bx
        fa = fb
        return
      end if
      ! v is beyond the limit, cutoff to the limit
    else if ((vx-cx)*(cx-bx) > Zero) then
      fv = bx+Lim*(cx-bx)
      Def = .false.
    end if
  end if
  ! unless the fit went beyond limits, use default step
  if (Def) then
    vx = cx+Ratio*(cx-bx)
    fv = f(q_a,q_b,dipole_a,dipole_b,r_a,r_b,vx,Absolute)
  end if
  ax = bx
  bx = cx
  cx = vx
  fa = fb
  fb = fc
  fc = fv
end do

end subroutine MnBrak2
