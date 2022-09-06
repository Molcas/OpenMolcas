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

subroutine Setup_Aux(nIrrep,nBas,nShell,nShell_Aux,nSO,TMax,CutOff,nij_Shell,nBas_Aux,nChV,iTOffs)

use iSD_data
use SOAO_Info, only: iSOInf
use j12, only: ShlSO, SOShl, nBasSh, iSSOff, iShij

implicit real*8(a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "setup.fh"
#include "nsd.fh"
integer nBas(0:nIrrep-1), nBas_Aux(0:nIrrep-1), nChV(0:nIrrep-1), iTOffs(3,0:nIrrep-1)
real*8 TMax(nShell,nShell)

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(6,*) 'Setup_Aux:nIrrep:   ',nIrrep
write(6,*) 'Setup_Aux:nBas:     ',nBas
write(6,*) 'Setup_Aux:nBas_Aux: ',nBas_Aux
write(6,*) 'Setup_Aux:nChV:     ',nChV
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
nSO = 0
nSO_Aux = 0
do iIrrep=0,nIrrep-1
  nSO = nSO+nBas(iIrrep)
  nSO_Aux = nSO_Aux+nBas_Aux(iIrrep)
end do

call mma_allocate(SOShl,nSO+nSO_Aux,Label='SOShl')
call mma_allocate(ShlSO,nSO+nSO_Aux,Label='ShlSO')
call mma_allocate(nBasSh,[0,nIrrep-1],[1,nShell+nShell_Aux],Label='nBasSh')
!                                                                      *
!***********************************************************************
!                                                                      *
do iSO=1,nSO+nSO_Aux
  iCnttp = iSOInf(1,iSO)
  iCnt = iSOInf(2,iSO)
  iAng = iSOInf(3,iSO)
  !write(6,*) 'iCnttp,iCnt,iAng=',iCnttp,iCnt,iAng

  ! Find the Shell from which this basis function is derived.

  do iSkal=1,nShell+nShell_Aux
    jCnttp = iSD(13,iSkal)
    jCnt = iSD(14,iSkal)
    jAng = iSD(1,iSkal)
    if ((jCnttp == iCnttp) .and. (jCnt == iCnt) .and. (jAng == iAng)) then
      SOShl(iSO) = iSkal
      !write(6,*) 'Found in shell=',iSkal
      Go To 99
    end if
  end do
99 continue
end do
!call iVcPrt('SOShl',' ',SOShl,nSO+nSO_Aux)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the number of effective shell pairs.

TMax_ij = Zero
do iSkal=1,nShell
  do jSkal=1,iSkal
    TMax_ij = max(TMax_ij,TMax(iSkal,jSkal))
  end do
end do

nij_Shell = 0
do iSkal=1,nShell
  do jSkal=1,iSkal
    if (TMax(iSkal,jSkal)*TMax_ij >= CutOff) then
      nij_Shell = nij_Shell+1
    end if
  end do
end do
call mma_allocate(iShij,2,nij_Shell,Label='iShij')

ij_Shell = 0
do iSkal=1,nShell
  do jSkal=1,iSkal
    if (TMax(iSkal,jSkal)*TMax_ij >= CutOff) then
      ij_Shell = ij_Shell+1
      iShij(1,ij_Shell) = iSkal
      iShij(2,ij_Shell) = jSkal
#     ifdef _DEBUGPRINT_
      write(6,*) 'ij_Shell,iSkal,jSkal=',ij_Shell,iSkal,jSkal
#     endif
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(iSSOff,[0,nIrrep-1],[0,nIrrep-1],[1,nij_Shell],Label='iSSOff')
!                                                                      *
!***********************************************************************
!                                                                      *
call Setup_Aux_Inner(SOShl,nSO+nSO_Aux,ShlSO,nBasSh,nShell+nShell_Aux,nIrrep,nBas,iSSOff,nij_Shell,iShij,nBas_Aux,nChV,iTOffs)
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(6,*) 'nij_Shell=',nij_Shell
write(6,*)
do ij_Shell=1,nij_Shell
  write(6,*) iShij(1,ij_Shell),iShij(2,ij_Shell)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Setup_Aux
