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

subroutine AddDPTC(nDPTC,nDSUM,DPTC,DSUM)

use caspt2_module, only: NSYM, NFRO, NORB, NBAS
use Constants, only: Half
use definitions, only: wp, iwp

implicit none

integer(kind=iwp), intent(in) :: nDPTC, nDSUM
real(kind=wp), intent(in) :: DPTC(nDPTC)
real(kind=wp), intent(inout) :: DSUM(nDSUM)

integer(kind=iwp) :: iMO1, iMO2, iSym, nOrbI1, nOrbI2, iOrb0, iOrb1, jOrb0, jOrb1, iOrb, jOrb
real(kind=wp) :: Val

iMO1 = 1
iMO2 = 1
do iSym=1,nSym
  nOrbI1 = nOrb(iSym)
  nOrbI2 = nBas(iSym)!-nDel(iSym)
  if (nOrbI2 > 0) then
    !! Add active orbital density
    !! Probably incorrect if symmetry
    do iOrb0=1,nOrb(iSym)
      iOrb1 = nFro(iSym)+iOrb0
      do jOrb0=1,nOrb(iSym)
        jOrb1 = nFro(iSym)+jOrb0
        DSUM(iMO1+iOrb0-1+nOrbI1*(jOrb0-1)) = DSUM(iMO1+iOrb0-1+nOrbI1*(jOrb0-1))+DPTC(iMO2+iOrb1-1+nOrbI2*(jOrb1-1))
      end do
    end do
    !! Symmetrize DSUM
    do iOrb=1,nOrb(iSym)
      do jOrb=1,iOrb-1
        Val = (DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1))+DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)))*Half
        DSUM(iMO1+iOrb-1+nOrbI1*(jOrb-1)) = Val
        DSUM(iMO1+jOrb-1+nOrbI1*(iOrb-1)) = Val
      end do
    end do
  end if
  iMO1 = iMO1+nOrbI1*nOrbI1
  iMO2 = iMO2+nOrbI2*nOrbI2
end do

end subroutine AddDPTC
