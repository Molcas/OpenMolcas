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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

subroutine WriteDip(DipGrad,Modes,Title,nOsc)

real*8 DipGrad(3,nOsc)
integer Modes(nOsc)
character Title*(*)
#include "inout.fh"

write(6,*)
write(6,*)
write(6,'(a2,a)') ' ',Title
write(6,'(a2,a)') ' ','=============================================='
write(6,'(a2,a)') ' ',' mode          X           Y           Z      '
write(6,'(a2,a)') ' ','----------------------------------------------'
do i=1,nOsc
  write(6,'(a3,i2,a1,a3,3f12.5)') ' ',Modes(i),'.',' ',(DipGrad(j,i),j=1,3)
end do
write(6,'(a2,a)') ' ','=============================================='
write(6,*)

end subroutine WriteDip
!####
subroutine IntCalcHeader()

#include "inout.fh"

write(6,*)
write(6,*)
write(6,'(a27,a)') ' ',' ================================================='
write(6,'(a27,a)') ' ','|                                                 |'
write(6,'(a27,a)') ' ','|             Intensity calculation               |'
write(6,'(a27,a)') ' ','|                                                 |'
write(6,'(a27,a)') ' ',' ================================================='
write(6,*)

end subroutine IntCalcHeader
!####
subroutine ExpPointHeader()

#include "inout.fh"

write(6,*)
write(6,*)
write(6,'(a27,a)') ' ',' ================================================='
write(6,'(a27,a)') ' ','|                                                 |'
write(6,'(a27,a)') ' ','|            Expansion point geometry             |'
write(6,'(a27,a)') ' ','|                                                 |'
write(6,'(a27,a)') ' ',' ================================================='
write(6,*)

end subroutine ExpPointHeader
!####
subroutine ISCHeader()

#include "inout.fh"

write(6,*)
write(6,*)
write(6,'(a27,a)') ' ',' ================================================='
write(6,'(a27,a)') ' ','|                                                 |'
write(6,'(a27,a)') ' ','|      InterSystem Crossing rate calculation      |'
write(6,'(a27,a)') ' ','|                                                 |'
write(6,'(a27,a)') ' ',' ================================================='
write(6,*)

end subroutine ISCHeader
!#####
subroutine WriteHeader(Title)
!  Purpose:
!    Write header and title to logfile.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

character*80 Title
#include "inout.fh"

write(6,*)

! Write title of project.
write(6,*)
write(6,*)
write(6,*)
write(6,'(A,A)') '  Title : ',Title
write(6,'(A)') '  -------'
write(6,*)

end subroutine WriteHeader

subroutine WrMold(FName,NumOfAt,AtomLbl,AtCoord,NumInt,HarmFreq,QMat)
! Open file named FName, and create a MOLDEN input file.

!implicit none
#include "Constants_mula.fh"
character*(*) FName
integer NumOfAt, NumInt
character*(*) AtomLbl(NumOfAt)
real*8 AtCoord(3,NumOfAt), HarmFreq(NumInt),QMat(3,NumOfAt,NumInt)
real*8 DMax, D2, DispMx, Factor
integer i, iInt, iAtom
character*2 AtName
#include "inout.fh"

call molcas_open(9,Fname)
!open(9,file=FName,status='UNKNOWN')
write(9,*) '[MOLDEN FORMAT]'

! Harmonic frequencies:
write(9,*) '[N_FREQ]'
write(9,*) NumInt
write(9,*) '[FREQ]'
do iInt=1,NumInt
  write(9,'(1X,F10.3)') HarToRcm*HarmFreq(iInt)
end do
! VV: Not implemented???
write(9,*) '[INT]'
do iInt=1,NumInt
  write(9,'(1X,F10.3)') 0.0
end do
! Atom coordinates:
write(9,*) '[NATOM]'
write(9,*) NumOfAt
write(9,*) '[FR-COORD]'
do iAtom=1,NumOfAt
  AtName = AtomLbl(iAtom)(1:2)
  if ((ichar('0') <= ichar(AtName(2:2))) .and. (ichar(AtName(2:2)) <= ichar('9'))) then
    AtName(2:2) = ' '
  end if
  write(9,'(1x,A2,3F16.8)') AtName,(AtCoord(i,iAtom),i=1,3)
end do

! Cartesian displacement coordinates:
! Scale such that max displacement is 0.2 A.U.
DMax = 0.2d0
write(9,*) '[FR-NORM-COORD]'
do iInt=1,NumInt
  write(9,*) 'vibration ',iInt
  DispMx = 0.0d0
  do iAtom=1,NumOfAt
    D2 = 0.0d0
    do i=1,3
      D2 = D2+QMat(i,iAtom,iInt)**2
    end do
    DispMx = max(DispMx,sqrt(D2))
  end do
  Factor = DMax/DispMx
  do iAtom=1,NumOfAt
    write(9,'(1X,3F16.8)') (Factor*QMat(i,iAtom,iInt),i=1,3)
  end do
end do

close(9)

end subroutine WrMold
