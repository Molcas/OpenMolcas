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

subroutine PtAdd(H0,Ovlp,RR,nSize,Temp,nTemp)
!***********************************************************************
!                                                                      *
!     Objective: Construct the modified Hamiltonian                    *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "input.fh"
real*8 H0(nSize), Ovlp(nSize), RR(nSize), Temp(nTemp)
character*8 Label
logical Debug
data Debug/.false./
dimension idum(1)

!----------------------------------------------------------------------*
!                                                                      *
!     Start procedure                                                  *
!     Read nuclear attraction and kinteic energy integrals.            *
!     Combine them to generate the one-electron Hamiltonian.           *
!     Finally read the overlap matrix.                                 *
!                                                                      *
!----------------------------------------------------------------------*

iOpt1 = 1
iOpt2 = 2
iComp = 1
iSyLbl = nSym
if (LCumulate) then
  Label = 'OneHam  '
  write(6,*)
  write(6,*) 'Adding perturbation cumulatively'
  write(6,*)
else
  Label = 'OneHam 0'
end if
iRc = -1
call iRdOne(iRc,iOpt1,Label,iComp,idum,iSyLbl)
nInts = idum(1)
if (iRc /= 0) then
  write(6,*) 'PtAdd: Error reading ONEINT'
  write(6,'(A,A)') 'Label=',Label
  call Abend()
end if
if (nInts+4 /= nSize) then
  write(6,*) 'PtAdd: nInts+4.ne.nSize',nInts+4,nSize
  call Abend()
end if
iRc = -1
call RdOne(iRc,iOpt2,Label,iComp,H0,iSyLbl)
if (Debug) then
  call PrDiOp('One Hamiltonian intgrl',nSym,nBas,H0)
  write(6,*) 'PotNuc=',H0(nInts+4)
end if

!----------------------------------------------------------------------*
!     Loop over all possible commands and branch to "special purpose"  *
!     subroutines to add perturbations.                                *
!----------------------------------------------------------------------*

call PtRela(H0,Ovlp,RR,nSize,Temp,nTemp)
call PtDipo(H0,Ovlp,RR,nSize,Temp,nTemp)
call PtQuad(H0,Ovlp,RR,nSize,Temp,nTemp)
call PtOkt0(H0,Ovlp,RR,nSize,Temp,nTemp)
call PtEfld(H0,Ovlp,RR,nSize,Temp,nTemp)
call PtEfgr(H0,Ovlp,RR,nSize,Temp,nTemp)
call PtGLbl(H0,Ovlp,RR,nSize,Temp,nTemp)

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
  write(6,*) 'PotNuc=',H0(nInts+4)
end if
iRc = -1
iOpt = 0
iComp = 1
Label = 'OneHam  '
call WrOne(iRc,iOpt,Label,iComp,H0,iSyLbl)
if (iRc /= 0) then
  write(6,*) 'PtAdd: Error writing to ONEINT'
  write(6,'(A,A)') 'Label=',Label
  call Abend()
end if
!call Put_PotNuc(H0(nInts+4))
call Put_dScalar('PotNuc',H0(nInts+4))

!----------------------------------------------------------------------*
!     Normal Exit                                                      *
!----------------------------------------------------------------------*

return

end subroutine PtAdd
