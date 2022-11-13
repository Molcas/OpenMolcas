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

subroutine PtAdd(H0,RR,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian                    *
!                                                                      *
!***********************************************************************

use FFPT_Global, only: LCumulate, nBas, nSym, ComStk
use OneDat, only: sNoOri, sOpSiz
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nSize, nTemp
real(kind=wp), intent(out) :: H0(nSize)
real(kind=wp), intent(inout) :: RR(nSize), Temp(nTemp)
character(len=8) :: Label
integer(kind=iwp) :: idum(1), iComp, iOpt, iOpt1, iOpt2, iRc, iSyLbl, nInts
logical(kind=iwp), parameter :: Debug = .false.

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Read nuclear attraction and kinteic energy integrals.            *
!     Combine them to generate the one-electron Hamiltonian.           *
!     Finally read the overlap matrix.                                 *
!                                                                      *
!----------------------------------------------------------------------*

iOpt1 = ibset(0,sOpSiz)
iOpt2 = ibset(0,sNoOri)
iComp = 1
iSyLbl = nSym
if (LCumulate) then
  Label = 'OneHam  '
  write(u6,*)
  write(u6,*) 'Adding perturbation cumulatively'
  write(u6,*)
else
  Label = 'OneHam 0'
end if
iRc = -1
call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
nInts = idum(1)
if (iRc /= 0) then
  write(u6,*) 'PtAdd: Error reading ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
if (nInts+4 /= nSize) then
  write(u6,*) 'PtAdd: nInts+4.ne.nSize',nInts+4,nSize
  call Abend()
end if
iRc = -1
call RdOne(iRc,iOpt2,Label,iComp,H0,iSyLbl)
if (Debug) then
  call PrDiOp('One Hamiltonian intgrl',nSym,nBas,H0)
  write(u6,*) 'PotNuc=',H0(nInts+4)
end if

!----------------------------------------------------------------------*
!     Loop over all possible commands and branch to "special purpose"  *
!     subroutines to add perturbations.                                *
!----------------------------------------------------------------------*

call PtRela(H0,nSize,Temp,nTemp)
call PtDipo(H0,nSize,Temp,nTemp)
call PtQuad(H0,RR,nSize,Temp,nTemp)
call PtOkt0(H0,RR,nSize,Temp,nTemp)
call PtEfld(H0,nSize,Temp,nTemp)
call PtEfgr(H0,nSize,Temp,nTemp)
call PtGLbl(H0,nSize,Temp,nTemp)

!----------------------------------------------------------------------*
!     If the user have requested a local (a la LoProp) perturbation    *
!     then make some modifications to the perturbation matrix.         *
!----------------------------------------------------------------------*

if (ComStk(4,0,0,0)) call SelectLoc(H0,nSize)

!----------------------------------------------------------------------*
!     Terminate procedure                                              *
!----------------------------------------------------------------------*

if (Debug) then
  call PrDiOp('Core Hamiltonian',nSym,nBas,H0)
  write(u6,*) 'PotNuc=',H0(nInts+4)
end if
iRc = -1
iOpt = 0
iComp = 1
Label = 'OneHam  '
call WrOne(iRc,iOpt,Label,iComp,H0,iSyLbl)
if (iRc /= 0) then
  write(u6,*) 'PtAdd: Error writing to ONEINT'
  write(u6,'(A,A)') 'Label=',Label
  call Abend()
end if
!call Put_PotNuc(H0(nInts+4))
call Put_dScalar('PotNuc',H0(nInts+4))

!----------------------------------------------------------------------*
!     Normal Exit                                                      *
!----------------------------------------------------------------------*

return

end subroutine PtAdd
