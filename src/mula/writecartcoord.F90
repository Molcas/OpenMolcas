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

character*4 AtomLbl(NumOfAt)
real*8 Coord(3,NumOfAt)
real*8 Mass(NumOfAt)
#include "inout.fh"

! Write labels, coordinates and masses to log file.
write(6,*)
write(6,*)
write(6,'(a1,a)') ' ','Cartesian coordinates (in bohr) and masses (in u)'
write(6,*) ('====',i=1,17)
write(6,'(a2,a)') ' ','Atom         x             y             z                Mass'
write(6,*) ('----',i=1,17)
do j=1,NumOfAt
  write(6,'(a2,a4,3f14.8,f20.8)') ' ',AtomLbl(j),(Coord(i,j),i=1,3),Mass(j)
end do
write(6,*) ('====',i=1,17)
write(6,*)
write(6,*)

end subroutine WriteCartCoord
!####
subroutine WriteIntCoord(InterVec,AtomLbl,xvec,NumInt)
!  Purpose:
!    Write internal coordinates to log file.

#include "Constants_mula.fh"
!integer InterVec(NumInt)
integer InterVec(*)
character*4 AtomLbl(NumInt)
character*128 Line
real*8 xvec(NumInt)
real*8 const
#include "inout.fh"

! Internal coordinates at equilibrium.
write(6,*)
write(6,*)
write(6,*)
write(6,'(a1,a)') ' ','Internal coordinates at equilibrium'
write(6,*) ('====',i=1,16)
write(6,'(a2,a)') ' ','Distances :                            bohr          aangstrom'
write(6,'(a2,a)') ' ','Angles    :                           radians         degrees'
write(6,*) ('----',i=1,16)
k = 1
do j=1,NumInt
  IntType = InterVec(k)
  if (IntType == 1) then
    ! Bond Stretching.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    write(Line,fmt='(A2,A,A,A,A,A)') ' ','Bond    ',AtomLbl(i1),'- ',AtomLbl(i2),'            '
    k = k+3
  end if
  if (IntType == 2) then
    ! Valence Angle Bending.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    write(Line,fmt='(A2,A,A,A,A,A,A,A)') ' ','Angle   ',AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),'      '
    k = k+4
  end if
  if (IntType == 3) then
    ! Linear Valence Angle.
    i1 = InterVec(k+1)
    i2 = InterVec(k+2)
    i3 = InterVec(k+3)
    write(Line,fmt='(A2,A,A,A,A,A,A,A)') ' ','LinAng  ',AtomLbl(i1),'- ',AtomLbl(i2),'- ',AtomLbl(i3),'      '
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
    const = 180.0d0/rpi
  end if
  write(6,'(A,A1,F15.8,F16.8)') Line(1:32),' ',xvec(j),xvec(j)*const
end do
write(6,*) ('====',i=1,16)
write(6,*)
write(6,*)

end subroutine WriteIntCoord
