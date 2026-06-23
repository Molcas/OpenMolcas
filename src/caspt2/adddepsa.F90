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

subroutine AddDEPSA(nDPT2,nAshT,DPT2,DEPSA)

use caspt2_module, only: NAES, NASH, NBAS, NDEL, NFRO, NISH, NORB, NSYM
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nDPT2, nAshT
real(kind=wp), intent(inout) :: DPT2(nDPT2)
real(kind=wp), intent(in) :: DEPSA(nAshT,nAshT)
integer(kind=iwp) :: iMO1, iMO2, iOrb, iOrb0, iOrb2, iSym, itabs, iuabs, jOrb, jOrb0, jOrb2, nOrbI1, nOrbI2
real(kind=wp) :: Val

iMO1 = 1
iMO2 = 1
do iSym=1,nSym
  nOrbI1 = nOrb(iSym)
  nOrbI2 = nBas(iSym)-nDel(iSym)
  if (nOrbI2 > 0) then
    !! Add active orbital density
    !! Probably incorrect if symmetry
    do iOrb0=1,nAsh(iSym)
      ! iOrb1 = nIsh(iSym)+iOrb0
      iOrb2 = nFro(iSym)+nIsh(iSym)+iOrb0
      itabs = iOrb0+NAES(iSym)
      do jOrb0=1,nAsh(iSym)
        ! jOrb1 = nIsh(iSym)+jOrb0
        jOrb2 = nFro(iSym)+nIsh(iSym)+jOrb0
        iuabs = jOrb0+NAES(iSym)
        DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1)) = DPT2(iMO2+iOrb2-1+nOrbI2*(jOrb2-1))+DEPSA(itabs,iuabs)
      end do
    end do
    !! Symmetrize DPT2 (for shift)
    do iOrb=1,nBas(iSym)-nDel(iSym)
      do jOrb=1,iOrb-1
        Val = (DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1))+DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)))*Half
        DPT2(iMO2+iOrb-1+nOrbI2*(jOrb-1)) = Val
        DPT2(iMO2+jOrb-1+nOrbI2*(iOrb-1)) = Val
      end do
    end do
  end if
  iMO1 = iMO1+nOrbI1*nOrbI1
  iMO2 = iMO2+nOrbI2*nOrbI2
end do

end subroutine AddDEPSA
