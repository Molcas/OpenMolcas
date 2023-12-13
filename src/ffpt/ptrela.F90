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

subroutine PtRela(H0,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Add the relativitic perturbation operator to          *
!                the one-electron Hamiltonian                          *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: nBas, nSym, ComStk, ComVal
use OneDat, only: sNoOri, sOpSiz
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize, nTemp
real(kind=wp), intent(inout) :: H0(nSize), Temp(nTemp)
character(len=8) :: Label
character(len=20) :: PriLbl
integer(kind=iwp) :: idum(1), iComp, iOpt1, iOpt2, iRc, iSyLbl, nInts
real(kind=wp) :: Alpha
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Load mass-velocity and one-electron Darwin contact term.         *
!     Add them to the one-electron Hamiltonian.                        *
!                                                                      *
!----------------------------------------------------------------------*

if (.not. ComStk(2,5,0,1)) then
  return
end if

Label = 'MassVel '
iRc = -1
iOpt1 = ibset(0,sOpSiz)
iOpt2 = ibset(0,sNoOri)
iSyLbl = 0
iComp = 1
Alpha = ComVal(2,5,0,1)
call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
nInts = idum(1)
if (iRc /= 0) call error()
call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
if (iRc /= 0) call error()
call daxpy_(nInts,Alpha,Temp,1,H0,1)
H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
if (Debug) then
  write(u6,'(6X,A,F8.6)') 'weight =',Alpha
  PriLbl = 'Mass-Velocity term  '
  write(PriLbl(19:20),'(I2)') iComp
  call PrDiOp(PriLbl,nSym,nBas,Temp)
end if
Label = 'Darwin  '
iRc = -1
iOpt1 = ibset(0,sOpSiz)
iOpt2 = ibset(0,sNoOri)
iSyLbl = 0
iComp = 1
call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
nInts = idum(1)
if (iRc /= 0) call error()
call RdOne(iRc,iOpt2,Label,iComp,Temp,iSyLbl)
if (iRc /= 0) call error()
call daxpy_(nInts,Alpha,Temp,1,H0,1)
H0(nInts+4) = H0(nInts+4)-Alpha*Temp(nInts+4)
if (Debug) then
  write(u6,'(6X,A,F8.6)') 'weight =',Alpha
  PriLbl = '1el. Darwin term    '
  write(PriLbl(19:20),'(I2)') iComp
  call PrDiOp(PriLbl,nSym,nBas,Temp)
end if

!----------------------------------------------------------------------*
!     Normal Exit                                                      *
!----------------------------------------------------------------------*

return

contains

!----------------------------------------------------------------------*
!     Error Exit                                                       *
!----------------------------------------------------------------------*
subroutine error()

  write(u6,*) 'PtRela: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()

end subroutine error

end subroutine PtRela
