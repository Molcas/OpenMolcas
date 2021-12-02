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

subroutine MnBrak(ax,bx,cx,fa,fb,fc,f,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New, &
                  Scratch_Org,iPrint_Errors)

use Constants, only: Zero, One, Five, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: ax, bx
real(kind=wp), intent(out) :: cx, fa, fb, fc
! External function f and its arguments
real(kind=wp), external :: f
integer(kind=iwp), intent(in) :: ij, l, nij, lMax, nElem, nAtoms, nPert, iPrint_Errors
real(kind=wp), intent(in) :: rMP(nij,0:nElem), xnrMP(nij,nElem), EC(3,nij), R_ij(3), C_o_C(3)
real(kind=wp), intent(out) :: xrMP(nij,nElem), xxrMP(nij,nElem), Scratch_New(nij*(2+lMax+1)), Scratch_Org(nij*(2+lMax+1))
real(kind=wp), intent(inout) :: AC(3,nij)
real(kind=wp) :: vx, fv, coefA, coefB
logical(kind=iwp) :: Def
real(kind=wp), parameter :: Lim = 100.0_wp, Ratio = Half*(One+sqrt(Five)), Thr = 1.0e-20_wp

fa = f(ax,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
fb = f(bx,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
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
fc = f(cx,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
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
      fv = f(vx,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
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
      fv = f(vx,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
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
    fv = f(vx,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
  end if
  ax = bx
  bx = cx
  cx = vx
  fa = fb
  fb = fc
  fc = fv
end do

end subroutine MnBrak
