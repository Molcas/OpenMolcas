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

subroutine OITD(rK,isym,D,Dtmp,act)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: G1t, ipCM, ipMat, nA, nDens
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(in) :: rK(*)
integer(kind=iwp), intent(in) :: iSym
real(kind=wp), intent(out) :: D(*)
real(kind=wp), intent(_OUT_) :: Dtmp(nDens)
logical(kind=iwp), intent(in) :: act
integer(kind=iwp) :: iB, iS, jB, jS

DTmp(:) = Zero

! Note: even with NAC we set the inactive block,
! because this is the SA density, not the transition density
do iS=1,nSym
  do iB=1,nIsh(iS)
    Dtmp(1+(ipCM(iS)+(iB-1)*nOrb(iS)+iB-1)-1) = Two
  end do
end do
if (act) then
  do iS=1,nSym
    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        Dtmp(1+(ipCM(iS)+iB+nIsh(iS)+(jB+nIsh(iS)-1)*nOrb(iS)-1)-1) = G1t((iTri((nA(iS)+iB),(nA(iS)+jB))))
      end do
    end do
  end do
end if

do iS=1,nsym
  jS = Mul(iS,isym)
  if (nOrb(iS)*nOrb(jS) >= 1) then
    call DGEMM_('N','T',nOrb(iS),nOrb(jS),nOrb(iS),One,Dtmp(1+ipCM(iS)-1),nOrb(iS),rK(ipMat(jS,iS)),nOrb(jS),Zero,D(ipMat(iS,jS)), &
                nOrb(iS))
    call DGEMM_('T','N',nOrb(iS),nOrb(jS),nOrb(jS),-One,rK(ipMat(jS,iS)),nOrb(jS),Dtmp(1+ipCM(jS)-1),nOrb(jS),One,D(ipMat(iS,jS)), &
                nOrb(iS))
  end if
end do

end subroutine OITD
