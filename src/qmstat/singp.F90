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

subroutine SingP(nCalls,iQ_Atoms,StoreCoo)

use qmstat_global, only: Cordst, DelFi, DelR, DelX, nAtom, nCent, nMacro, nMicro, nPart, Qmeq
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: nCalls
integer(kind=iwp), intent(in) :: iQ_Atoms
real(kind=wp), intent(inout) :: StoreCoo(3,nCent,nPart)
integer(kind=iwp) :: iCent, Initial1, iPart, kaunter, nAllQm

if (nCalls == 0) then
  ! If this is first call, issue a warning.

  write(u6,*)
  write(u6,*)
  write(u6,*) '---->>>  WARNING  <<<----'
  write(u6,*)
  write(u6,*) 'You have specified that a set of single-point calculations are to be performed.'
  write(u6,*) 'This means that the input will be given to some extent a new meaning.'
  write(u6,*)

  ! Put coordinates in a new vector if first call.

  kaunter = 0
  do iPart=1,nPart
    do iCent=1,nCent
      kaunter = kaunter+1
      StoreCoo(:,iCent,iPart) = Cordst(:,kaunter)
    end do
  end do

  ! Put dummies that will be substituted for the qm-region.

  nAllQm = (((iQ_Atoms-1)/nAtom)+1)*nCent
  Cordst(:,1:nAllQm) = Zero

  ! Put the coordinates of first iteration.

  Cordst(:,nAllQm+1:nAllQm+nCent) = StoreCoo(:,:,1)

  ! Set new value on some variables.

  nMicro = 1
  nMacro = 1
  DelX = Zero
  DelFi = Zero
  DelR = Zero
  QmEq = .true.
  nPart = (nAllQm/nCent)+1
  write(u6,*)
  write(u6,*) 'Resetting for FIT:'
  write(u6,*) 'Number of macrosteps:',nMacro
  write(u6,*) 'Number of microsteps:',nMicro
  write(u6,*) 'No translation, rotation or radie modification.'
  write(u6,*) 'Take the QmEq path.'

else
  ! If not first call, then collect relevant coordinates.

  Initial1 = (((iQ_Atoms-1)/nAtom)+1)*nCent
  Cordst(:,Initial1+1:Initial1+nCent) = StoreCoo(:,:,nCalls+1)
end if

! Update nCalls.

nCalls = nCalls+1

return

end subroutine SingP
