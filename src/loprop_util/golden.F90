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

function Golden(ax,bx,cx,f,tol_x,tol_f,xmin,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New, &
                Scratch_Org,iPrint_Errors)

use Constants, only: One, Three, Five, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: Golden
real(kind=wp), intent(in) :: ax, bx, cx, tol_x, tol_f
real(kind=wp), intent(out) :: xmin
! External function f and its arguments
real(kind=wp), external :: f
integer(kind=iwp), intent(in) :: ij, l, nij, lMax, nElem, nAtoms, nPert, iPrint_Errors
real(kind=wp), intent(in) :: rMP(nij,0:nElem), xnrMP(nij,nElem), EC(3,nij), R_ij(3), C_o_C(3)
real(kind=wp), intent(out) :: xrMP(nij,nElem), xxrMP(nij,nElem), Scratch_New(nij*(2+lMax+1)), Scratch_Org(nij*(2+lMax+1))
real(kind=wp), intent(inout) :: AC(3,nij)
real(kind=wp) :: f2, f3, x1, x2, x3, x4
real(kind=wp), parameter :: Ratio = Half*(Three-sqrt(Five)), RM = One-Ratio

x1 = ax
x4 = cx
if (abs(cx-bx) > abs(bx-ax)) then
  x2 = bx
  x3 = RM*bx+Ratio*cx
else
  x2 = RM*bx+Ratio*ax
  x3 = bx
end if
f2 = f(x2,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
f3 = f(x3,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
do while ((abs(x4-x1) > tol_x*(abs(x1)+abs(x2))) .and. (abs(f3-f2) > tol_f*(abs(f2)+abs(f3))))
  if (f2 < f3) then
    x4 = x3
    x3 = x2
    f3 = f2
    x2 = RM*x3+Ratio*x1
    f2 = f(x2,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
  else
    x1 = x2
    x2 = x3
    f2 = f3
    x3 = RM*x2+Ratio*x4
    f3 = f(x3,rMP,xrMP,xxrMP,xnrMP,EC,AC,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org,iPrint_Errors)
  end if
end do
if (f2 < f3) then
  xmin = x2
  Golden = f2
else
  xmin = x3
  Golden = f3
end if

end function Golden
