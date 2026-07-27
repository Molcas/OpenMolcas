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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

! MO->AO or AO->MO transformation of 1-RDM
subroutine OLagTrf(mode,iSym,NBSQT,CMO,DPT2,DPT2AO,WRK)

use caspt2_module, only: NBAS, NDEL, NFRO, NORB
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: mode, iSym, NBSQT
real(kind=wp), intent(in) :: CMO(NBSQT)
real(kind=wp), intent(inout) :: DPT2(NBSQT), DPT2AO(NBSQT)
real(kind=wp), intent(out) :: WRK(NBSQT)
real(kind=wp) :: Val
integer(kind=iwp) :: iBas, iCMO, iMO, jBas, nBasI, nOrbI

!! Mode = 1: MO -> AO transformation
!! Mode = 2: AO -> MO transformation
iCMO = 1+sum(nBas(1:iSym-1)**2)
iMO = 1+sum((nOrb(1:iSym-1)+nFro(1:iSym-1))**2)

if (nOrb(iSym)+nFro(iSym) > 0) then
  nBasI = nBas(iSym)
  nOrbI = nBas(iSym)-nDel(iSym)
  if (Mode == 1) then
    !! MO -> AO
    call DGEMM_('N','N',nBasI,nOrbI,nOrbI,One,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,Zero,WRK,nBasI)
    call DGEMM_('N','T',nBasI,nBasI,nOrbI,One,WRK,nBasI,CMO(iCMO),nBasI,Zero,DPT2AO(iCMO),nBasI)
    !! Symmetrize, just in case
    do iBas=1,nBasI
      do jBas=1,iBas-1
        Val = (DPT2AO(iCMO+iBas-1+nBasI*(jBas-1))+DPT2AO(iCMO+jBas-1+nBasI*(iBas-1)))*Half
        DPT2AO(iCMO+iBas-1+nBasI*(jBas-1)) = Val
        DPT2AO(iCMO+jBas-1+nBasI*(iBas-1)) = Val
      end do
    end do
  else if (Mode == 2) then
    !! AO -> MO
    call DGEMM_('T','N',nOrbI,nBasI,nBasI,One,CMO(iCMO),nBasI,DPT2AO(iCMO),nBasI,Zero,WRK,nOrbI)
    call DGEMM_('N','N',nOrbI,nOrbI,nBasI,One,WRK,nOrbI,CMO(iCMO),nBasI,Zero,DPT2(iMO),nOrbI)
  end if
end if

return

end subroutine OLagTrf
