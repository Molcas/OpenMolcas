!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine abc_axes(cryst,coord,xyz,abc,Do_option,iReturn)
! this Subroutine performs a transformation of the main axes
! (magnetic, anisotropic etc.) from a xyz system in
! crystallographic "abc" system, and vice-versa;
!    Do_option = 1 =>  transform from xyz to abc
!    Do_option = 2 =>  transform from abc to xyz
!    If Do_option has other value, abort
!    coord(3)-- Cartesian coordinates of the main magnetic center

use Constants, only: Zero, One, Two, deg2rad
use Definitions, only: u6

implicit none
integer, intent(in) :: Do_option
real(kind=8), intent(inout) :: xyz(3,3), abc(3,3)
real(kind=8), intent(in) :: cryst(6), coord(3)
integer, intent(out) :: iReturn
! local variables:
integer :: i
real(kind=8) :: a, b, c, al, bt, gm, cal, cbt, cgm, sgm, v, x, y, z
real(kind=8) :: xyz2(3,3), pX(3), pY(3), pZ(3)

! initializations
a = cryst(1)
b = cryst(2)
c = cryst(3)
al = cryst(4)*deg2rad
bt = cryst(5)*deg2rad
gm = cryst(6)*deg2rad
cal = cos(al)
cbt = cos(bt)
cgm = cos(gm)
sgm = sin(gm)

X = Zero
Y = Zero
Z = Zero
pX(:) = Zero
pY(:) = Zero
pZ(:) = Zero
xyz2(:,:) = Zero

v = sqrt(One-cal*cal-cbt*cbt-cgm*cgm+Two*cal*cbt*cgm)

if (Do_option == 1) then
  abc = Zero
  do i=1,3
    xyz2(1,i) = xyz(1,i)+coord(1)
    xyz2(2,i) = xyz(2,i)+coord(2)
    xyz2(3,i) = xyz(3,i)+coord(3)
  end do

  do i=1,3
    X = xyz2(1,i)
    Y = xyz2(2,i)
    Z = xyz2(3,i)

    pX(1) = One/a
    pY(1) = -cgm/(a*sgm)
    pZ(1) = ((cal*cgm-cbt)/(a*v*sgm))

    pX(2) = Zero
    pY(2) = One/(b*sgm)
    pZ(2) = (cbt*cgm-cal)/(b*v*sgm)

    pX(3) = Zero
    pY(3) = Zero
    pZ(3) = sgm/(c*v)

    abc(1,i) = pX(1)*X+pY(1)*Y+pZ(1)*Z
    abc(2,i) = pX(2)*X+pY(2)*Y+pZ(2)*Z
    abc(3,i) = pX(3)*X+pY(3)*Y+pZ(3)*Z
  end do

else if (Do_option == 2) then

  xyz = Zero
  do i=1,3
    X = abc(1,i)*a
    Y = abc(2,i)*b
    Z = abc(3,i)*c

    pX(1) = One
    pY(1) = cgm
    pZ(1) = cbt

    pX(2) = Zero
    pY(2) = sgm
    pZ(2) = (cal-cbt*cgm)/sgm

    pX(3) = Zero
    pY(3) = Zero
    pZ(3) = v/sgm

    xyz(1,i) = X+Y*cgm+Z*cbt
    xyz(2,i) = Y*sgm+Z*((cal-cbt*cgm)/sgm)
    xyz(3,i) = Z*v/sgm

    xyz(1,i) = pX(1)*X+pY(1)*Y+pZ(1)*Z
    xyz(2,i) = pX(2)*X+pY(2)*Y+pZ(2)*Z
    xyz(3,i) = pX(3)*X+pY(3)*Y+pZ(3)*Z
  end do
else
  write(u6,'(A)') 'the Do_option is not specified. '
  write(u6,'(A)') 'the program continues without ABCC option'
  iReturn = 1
end if

return

end subroutine abc_axes
