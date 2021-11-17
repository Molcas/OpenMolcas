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
write(u6,*) ('====',i=1,17)
write(u6,'(a2,a)') ' ','Atom         x             y             z                Mass'
write(u6,*) ('----',i=1,17)
do i=1,NumOfAt
  write(u6,'(a2,a4,3f14.8,f20.8)') ' ',AtomLbl(i),Coord(:,i),Mass(i)
end do
write(u6,*) ('====',i=1,17)
write(u6,*)
write(u6,*)

end subroutine WriteCartCoord
!####
subroutine WriteIntCoord(InterVec,AtomLbl,xvec,NumInt)
!  Purpose:
!    Write internal coordinates to log file.

use Constants, only: One, Angstrom, deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: InterVec(*), NumInt
character(len=4), intent(in) :: AtomLbl(NumInt)
real(kind=wp), intent(in) :: xvec(NumInt)
integer(kind=iwp) :: i, i1, i2, i3, i4, IntType, j, k
real(kind=wp) :: const
character(len=128) :: Line

! Internal coordinates at equilibrium.
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,'(a1,a)') ' ','Internal coordinates at equilibrium'
write(u6,*) ('====',i=1,16)
write(u6,'(a2,a)') ' ','Distances :                            bohr           angstrom'
write(u6,'(a2,a)') ' ','Angles    :                           radians         degrees'
write(u6,*) ('----',i=1,16)
k = 1
do j=1,NumInt
  IntType = InterVec(k)
  if (IntType == 1) then
    ! Bond Stretching.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    write(Line,fmt='(A2,A,A,A,A)') ' ','Bond    ',AtomLbl(i1),'- ',AtomLbl(i2)
    k = k+3
  end if
  if (IntType == 2) then
    ! Valence Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    write(Line,fmt='(A2,A,A,A,A,A,A)') ' ','Angle   ',AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3)
    k = k+4
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    write(Line,fmt='(A2,A,A,A,A,A,A)') ' ','LinAng  ',AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3)
    k = k+4
  end if
  if (IntType == 4) then
    ! Torsion.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    write(Line,fmt='(A2,A,A,A,A,A,A,A,A)') ' ','Torsion ',AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),'- ',AtomLbl(i4)
    k = k+5
  end if
  if (IntType == 5) then
    ! Out of Plane Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    i4 = InterVec(k+4)
    write(Line,fmt='(A2,A,A,A,A,A,A,A,A)') ' ','OutOfPl ',AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),'- ',AtomLbl(i4)
    k = k+5
  end if
  if (intType == 1) then
    const = Angstrom
  else
    const = One/deg2rad
  end if
  write(u6,'(A,A1,F15.8,F16.8)') Line(1:32),' ',xvec(j),xvec(j)*const
end do
write(u6,*) ('====',i=1,16)
write(u6,*)
write(u6,*)

end subroutine WriteIntCoord
