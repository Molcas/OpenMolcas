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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine DPT2_Trf(NBSQT,nAshT,DPT2,DPT2AO,CMO,DEPSA,DSUM)

use general_data, only: NASH
use caspt2_module, only: NBAS, NDEL, NFRO, NISH, NORB, NSYM
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT, nAshT
real(kind=wp), intent(inout) :: DPT2(NBSQT), DPT2AO(NBSQT), DSUM(NBSQT)
real(kind=wp), intent(in) :: CMO(NBSQT), DEPSA(nAshT,nAshT)
integer(kind=iwp) :: iAO, iCMO, iMO, iOrb, iOrb0, iSym, jOrb, jOrb0, nBasI, nOrbI
real(kind=wp) :: Val
real(kind=wp), allocatable :: WRK(:)

!! DPT2 transformation
!! Just transform DPT2 (in MO, block-squared) to DPT2AO (in AO,
!! block-squared). Also, for DPT2C which couples with the inactive
!! density matrix.
call mma_allocate(WRK,NBSQT,Label='WRK')

!! MO -> AO back transformation
iCMO = 1
iAO = 1
iMO = 1
do iSym=1,nSym
  iCMO = iCMO+nBas(iSym)*nFro(iSym)
  if (nORB(ISYM) > 0) then
    nBasI = nBas(iSym)
    nOrbI = nOrb(iSym)
    !! Add active orbital density
    do iOrb0=1,nAsh(iSym)
      iOrb = nIsh(iSym)+iOrb0
      ! iOrb2= nFro(iSym)+nIsh(iSym)+iOrb0
      do jOrb0=1,nAsh(iSym)
        jOrb = nIsh(iSym)+jOrb0
        ! jOrb2= nFro(iSym)+nIsh(iSym)+jOrb0
        DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))+DEPSA(iOrb0,jOrb0)
        DSUM(iMO+iOrb-1+nOrbI*(jOrb-1)) = DSUM(iMO+iOrb-1+nOrbI*(jOrb-1))+DEPSA(iOrb0,jOrb0)
      end do
    end do
    !! Symmetrize DPT2 (for shift)
    do iOrb=1,nOrb(iSym)
      do jOrb=1,iOrb
        Val = (DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))+DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)))*Half
        DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Val
        DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Val
      end do
    end do
    !! First, DPT2 -> DPT2AO
    call DGEMM_('N','N',nBasI,nOrbI,nOrbI,One,CMO(iCMO),nBasI,DPT2(iMO),nOrbI,Zero,WRK,nBasI)
    call DGEMM_('N','T',nBasI,nBasI,nOrbI,One,WRK,nBasI,CMO(iCMO),nBasI,Zero,DPT2AO(iAO),nBasI)
  end if
  iCMO = iCMO+nBas(iSym)*(nOrb(iSym)+nDel(iSym))
  iAO = iAO+nBasI*nBasI
  iMO = iMO+nBasI*nBasI
end do

call mma_deallocate(WRK)

end subroutine DPT2_Trf
