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

subroutine SingP(nCalls,iQ_Atoms,ipStoreCoo,nPart2)

implicit real*8(a-h,o-z)
#include "maxi.fh"
#include "qminp.fh"
#include "WrkSpc.fh"

if (nCalls == 0) then
  ! If this is first call, issue a warning.

  write(6,*)
  write(6,*)
  write(6,*) '---->>>  WARNING  <<<----'
  write(6,*)
  write(6,*) 'You have specified that a set of single-point calculations are to be preformed.'
  write(6,*) 'This means that the input will be given to some extent a new meaning.'
  write(6,*)

  ! Put coordinates in a new vector if first call.

  kaunter = 0
  nPart2 = nPart
  call GetMem('Store','Allo','Real',ipStoreCoo,nPart2*nCent*3)
  do iPart=1,nPart2
    do iCent=1,nCent
      kaunter = kaunter+1
      Work(ipStoreCoo+3*(kaunter-1)) = Cordst(kaunter,1)
      Work(ipStoreCoo+3*(kaunter-1)+1) = Cordst(kaunter,2)
      Work(ipStoreCoo+3*(kaunter-1)+2) = Cordst(kaunter,3)
    end do
  end do

  ! Put dummies that will be substituted for the qm-region.

  nAllQm = (((iQ_Atoms-1)/nAtom)+1)*nCent
  do i=1,nAllQm
    do j=1,3
      Cordst(i,j) = 0.0d0
    end do
  end do

  ! Put the coordinates of first iteration.

  do iCent=1,nCent
    Cordst(nAllQm+iCent,1) = Work(ipStoreCoo+3*(iCent-1))
    Cordst(nAllQm+iCent,2) = Work(ipStoreCoo+3*(iCent-1)+1)
    Cordst(nAllQm+iCent,3) = Work(ipStoreCoo+3*(iCent-1)+2)
  end do

  ! Set new value on some variables.

  nMicro = 1
  nMacro = 1
  DelX = 0
  DelFi = 0
  DelR = 0
  QmEq = .true.
  nPart = (nAllQm/nCent)+1
  write(6,*)
  write(6,*) 'Resetting for FIT:'
  write(6,*) 'Number of macrosteps:',nMacro
  write(6,*) 'Number of microsteps:',nMicro
  write(6,*) 'No translation, rotation or radie modification.'
  write(6,*) 'Take the QmEq path.'


else
  ! If not first call, then collect relevant coordinates.

  Initial1 = (((iQ_Atoms-1)/nAtom)+1)*nCent
  Initial2 = 3*nCent*nCalls-1
  do iCent=1,nCent
    Cordst(Initial1+iCent,1) = Work(ipStoreCoo+Initial2+(iCent-1)*3+1)
    Cordst(Initial1+iCent,2) = Work(ipStoreCoo+Initial2+(iCent-1)*3+2)
    Cordst(Initial1+iCent,3) = Work(ipStoreCoo+Initial2+(iCent-1)*3+3)
  end do
end if

! Update nCalls.

nCalls = nCalls+1

return

end subroutine SingP
