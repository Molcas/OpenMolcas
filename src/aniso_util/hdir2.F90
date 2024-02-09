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

subroutine hdir2(nP,L,dX,dY,dZ,Ang,iprint)
! angstep == steps in the angular distribution. It defines the number of points
!            in which M will be computed;
!      nP == number of points, ( nP = 360/angstep )
!       L == cartesian component of the magnetisation torque
!            If L=1 (i.e.X), rotation of the M occurs in the YZ plane
!            If L=2 (i.e.Y), rotation of the M occurs in the XZ plane
!            If L=3 (i.e.Z), rotation of the M occurs in the XY plane

use Constants, only: Zero, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nP, L, iprint
real(kind=wp), intent(out) :: dX(nP), dY(nP), dZ(nP), Ang(nP)
integer(kind=iwp) :: i
real(kind=wp) :: AngRad, AngStep

AngStep = 360.0_wp/real(nP-1,kind=wp)

if (L == 1) then
  dX(:) = Zero
  do i=1,nP
    AngRad = real(i-1,kind=wp)*AngStep*deg2rad
    Ang(i) = real(i-1,kind=wp)*AngStep
    dY(i) = cos(AngRad)
    dZ(i) = sin(AngRad)
  end do
else if (L == 2) then
  dY(:) = Zero
  do i=1,nP
    AngRad = real(i-1,kind=wp)*AngStep*deg2rad+122.625_wp*deg2rad
    Ang(i) = real(i-1,kind=wp)*AngStep
    dX(i) = cos(AngRad)
    dZ(i) = sin(AngRad)
  end do
else if (L == 3) then
  dZ(:) = Zero
  do i=1,nP
    AngRad = real(i-1,kind=wp)*AngStep*deg2rad
    Ang(i) = real(i-1,kind=wp)*AngStep
    dX(i) = cos(AngRad)
    dY(i) = sin(AngRad)
  end do
else
  write(u6,'(A   )') 'Error. Parametr L can take only Integer values 1, 2 or 3.'
  write(u6,'(A,I5)') 'Current value: L = ',L
end if

if (iprint > 2) then
  write(u6,'(A,I5)') 'Angular grid for Magnetization Torque, Cartesian Component =',L
  write(u6,'(2x,A,4x,A,5x,3(10X,A,10x))') 'Nr.','Angle','X','Y','Z'
  do i=1,nP
    write(u6,'(I4,F10.3,3x,3F21.14)') i,Ang(i),dX(i),dY(i),dZ(i)
  end do
end if

return

end subroutine hdir2
