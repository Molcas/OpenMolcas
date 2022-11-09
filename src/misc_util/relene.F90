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

subroutine RelEne(ErelMV,ErelDC,nSym,nBas,CMO,OCC,D,OP)
!***********************************************************************
!                                                                      *
!     Purpose: Compute relativistic correction to energy for a given   *
!              set of natural orbitals and occupation numbers          *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri, sOpSiz
use Constants, only: Zero, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: ErelMV, ErelDC
integer(kind=iwp), intent(in) :: nSym, nBas(*)
real(kind=wp), intent(in) :: CMO(*), OCC(*)
real(kind=wp), intent(_OUT_) :: D(*), OP(*)
integer(kind=iwp) :: iBas, iComp, ij, iOff1, iOff2, iOff3, iOpt, iOrb, iRc, iSyLbl, iSym, jBas, lOp, nBs
character(len=8) :: Label
real(kind=wp), external :: DDOT_

!----------------------------------------------------------------------*
! Compute 1-particle density matrix in AO basis                        *
!----------------------------------------------------------------------*
ij = 0
iOff1 = 0
iOff2 = 0
do iSym=1,nSym
  nBs = nBas(iSym)
  if (nBs == 0) cycle
  do iBas=1,nBs
    do jBas=1,iBas
      ij = ij+1
      D(ij) = Zero
      do iOrb=1,nBs
        iOff3 = iOff1+(iOrb-1)*nBs
        D(ij) = D(ij)+OCC(iOff2+iOrb)*CMO(iBas+iOff3)*CMO(jBas+iOff3)
      end do
      if (iBas /= jBas) D(ij) = Two*D(ij)
    end do
  end do
  iOff1 = iOff1+nBs*nBs
  iOff2 = iOff2+nBs
end do
!----------------------------------------------------------------------*
! Compute energy contributions                                         *
!----------------------------------------------------------------------*
lOp = 0
do iSym=1,nSym
  nBs = nBas(iSym)
  lOp = lOp+nTri_Elem(nBs)
end do
ErelMV = Zero
iRc = -1
iOpt = ibset(0,sOpSiz)
iComp = 1
Label = 'MassVel'
call RdOne(iRc,iOpt,Label,iComp,OP,iSyLbl)
if (iRc == 0) then
  iRc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  call RdOne(iRc,iOpt,Label,iComp,OP,iSyLbl)
  ErelMV = DDOT_(lOP,D,1,OP,1)
end if
ErelDC = Zero
iRc = -1
iOpt = ibset(0,sOpSiz)
iComp = 1
Label = 'Darwin'
call RdOne(iRc,iOpt,Label,iComp,OP,iSyLbl)
if (iRc == 0) then
  iRc = -1
  iOpt = ibset(ibset(0,sNoOri),sNoNuc)
  iComp = 1
  call RdOne(iRc,iOpt,Label,iComp,OP,iSyLbl)
  ErelDC = DDOT_(lOP,D,1,OP,1)
end if

!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine RelEne
