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

subroutine caspt2_grad_invaria2(NDPT2,nOLag,DPT2,OLag)
! Compute pseudo-density in the inactive and secondary orbital blocks for MRPT2 methods that are not invariant with respect to
! orbital rotations among inactive and secondary orbitals using the canonical condition of MOs
! See the IPEA-shift implementation for the active block non-invariance, and the sigma^P implementation for the non-invariance
! with respect to rotations in the internally contracted basis

use Constants, only: Half
use definitions, only: iwp, wp
use caspt2_module, only: nSym, nDel, nIsh, EPSI, nSsh, nAsh, EPSE, nBas, nFro

implicit none
integer(kind=iwp), intent(in) :: NDPT2, nOLag
real(kind=wp), intent(inout) :: DPT2(NDPT2)
real(kind=wp), intent(in) :: OLag(nOLag)
integer(kind=iwp) :: iMO, iSym, nOrbI, nFroI, nIshI, nAshI, nSshI, iOrb, jOrb
real(kind=wp) :: Tmp

iMO = 1
do iSym=1,nSym
  nOrbI = nBas(iSym)-nDel(iSym)
  nFroI = nFro(iSym)
  nIshI = nIsh(iSym)
  if (nIshI > 0) then
    do iOrb=nFroI+1,nFroI+nIshI
      do jOrb=iOrb+1,nFroI+nIshI
        Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))-OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))/(EPSI(iOrb-nFroI)-EPSI(jOrb-nFroI))
        DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Tmp
        DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Tmp
      end do
    end do
  end if
  nSshI = nSsh(iSym)
  nAshI = nAsh(iSym)
  if (nSshI > 0) then
    do iOrb=nOrbI-nSshI+1,nOrbI
      do jOrb=iOrb+1,nOrbI
        Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))-OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))/ &
              (EPSE(iOrb-nFroI-nIshI-nAshI)-EPSE(jOrb-nFroI-nIshI-nAshI))
        DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = Tmp
        DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = Tmp
      end do
    end do
  end if
  iMO = iMO+nOrbI*nOrbI
end do

end subroutine caspt2_grad_invaria2
