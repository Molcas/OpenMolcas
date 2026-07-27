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

subroutine OLagFro1(NBSQT,nOLag,DPT2,OLag)

use caspt2_global, only: FIFA_all
use caspt2_module, only: NBAS, NDEL, NFRO, NISH, NSYM
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NBSQT, nOLag
real(kind=wp), intent(inout) :: DPT2(NBSQT), OLag(nOLag)
integer(kind=iwp) :: iMO, iOrb, iSym, jOrb, nBasI, nFroI, nIshI, nOrbI
real(kind=wp) :: Tmp

iMO = 1
do iSym=1,nSym
  nOrbI = nBas(iSym)-nDel(iSym)
  nFroI = nFro(iSym)
  if ((nOrbI > 0) .and. (nFroI > 0)) then
    nIshI = nIsh(iSym)
    nBasI = nBas(iSym)
    !! Make sure that the frozen orbital derivative is zero
    !! (it does not appear in the PT2 energy)
    OLag(1:nOrbI*nFroI) = Zero
    do iOrb=1,nFroI
      do jOrb=nFroI+1,nFroI+nIshI
        Tmp = -Half*(OLag(iMO+iOrb-1+nOrbI*(jOrb-1))-OLag(iMO+jOrb-1+nOrbI*(iOrb-1)))/ &
              (FIFA_all(iOrb+nBasI*(iOrb-1))-FIFA_all(jOrb+nBasI*(jOrb-1)))
        DPT2(iMO+iOrb-1+nOrbI*(jOrb-1)) = DPT2(iMO+iOrb-1+nOrbI*(jOrb-1))+Tmp
        DPT2(iMO+jOrb-1+nOrbI*(iOrb-1)) = DPT2(iMO+jOrb-1+nOrbI*(iOrb-1))+Tmp
      end do
    end do
  end if
  iMO = iMO+nOrbI*nOrbI
end do
!write(u6,*) 'DPT2 after frozen orbital'
!call sqprt(dpt2,nbast)

end subroutine OLagFro1
