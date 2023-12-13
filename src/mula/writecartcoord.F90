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

subroutine WriteCartCoord(AtomLbl,Coord,Mass,NumOfAt)
!  Purpose:
!    Write cartesian coordinates to log file.

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NumOfAt
character(len=4), intent(in) :: AtomLbl(NumOfAt)
real(kind=wp), intent(in) :: Coord(3,NumOfAt), Mass(NumOfAt)
integer(kind=iwp) :: i

! Write labels, coordinates and masses to log file.
write(u6,*)
write(u6,*)
write(u6,'(a1,a)') ' ','Cartesian coordinates (in bohr) and masses (in u)'
write(u6,*) repeat('=',51)
write(u6,'(a2,a)') ' ','Atom         x             y             z                Mass'
write(u6,*) repeat('-',51)
do i=1,NumOfAt
  write(u6,'(a2,a4,3f14.8,f20.8)') ' ',AtomLbl(i),Coord(:,i),Mass(i)
end do
write(u6,*) repeat('=',51)
write(u6,*)
write(u6,*)

end subroutine WriteCartCoord
