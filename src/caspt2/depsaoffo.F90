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

subroutine DEPSAOffO(nOLag,nAshT,NBSQT,OLag,DEPSA,FIFA)

use general_data, only: NASH
use caspt2_module, only: NBAS, NDEL, NFRO, NISH, NSYM
use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nOLag, nAshT, NBSQT
real(kind=wp), intent(in) :: OLag(nOLag), FIFA(NBSQT)
real(kind=wp), intent(inout) :: DEPSA(nAshT,nAshT)
integer(kind=iwp) :: iAsh, iMO, iOrb, iSym, jAsh, jOrb, nAshI, nBasI, nFroI, nIshI, nOrbI
real(kind=wp) :: EigI, EigJ, OLagIJ, OLagJI, Tmp

! This is much easier; similar to the frozen core approximation.
! Corresponds to the first term in Eq. (70)

iMO = 1
do iSym=1,nSym
  nAshI = nAsh(iSym)
  nOrbI = 0
  if (nAshI /= 0) then
    nOrbI = nBas(iSym)-nDel(iSym)
    nFroI = nFro(iSym)
    nIshI = nIsh(iSym)
    nBasI = nBas(iSym)

    do iAsh=1,nAshI
      iOrb = iAsh+nFroI+nIshI
      EigI = FIFA(iMO+iOrb-1+nBasI*(iOrb-1))
      do jAsh=1,iAsh-1
        jOrb = jAsh+nFroI+nIshI
        EigJ = FIFA(iMO+jOrb-1+nBasI*(jOrb-1))
        OLagIJ = OLag(iMO+iOrb-1+nOrbI*(jOrb-1))
        OLagJI = OLag(iMO+jOrb-1+nOrbI*(iOrb-1))
        Tmp = -(OLagIJ-OLagJI)/(EigI-EigJ)*Half

        DEPSA(iAsh,jAsh) = DEPSA(iAsh,jAsh)+Tmp
        DEPSA(jAsh,iAsh) = DEPSA(jAsh,iAsh)+Tmp
      end do
    end do
  end if
  iMO = iMO+nOrbI*nOrbI
end do

end subroutine DEPSAOffO
