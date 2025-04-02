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
use MCLR_Data, only: G1t
use MCLR_Data, only: ipCM, ipMat, nA, nDens2
use input_mclr, only: nSym, nAsh, nIsh, nOrb
use Constants, only: Zero, One, Two

implicit none
integer iSym
real*8 rK(*), D(*), Dtmp(nDens2)
logical act
integer iS, iB, jB, jS

DTmp(:) = Zero

! Note: even with NAC we set the inactive block,
! because this is the SA density, not the transition density
do iS=1,nSym
  do iB=1,nIsh(iS)
    Dtmp(1+(ipCM(iS)+(ib-1)*nOrb(iS)+ib-1)-1) = Two
  end do
end do
if (act) then
  do iS=1,nSym
    do iB=1,nAsh(iS)
      do jB=1,nAsh(iS)
        Dtmp(1+(ipCM(iS)+ib+nIsh(is)+(jB+nIsh(is)-1)*nOrb(is)-1)-1) = G1t((iTri((nA(is)+ib),(nA(is)+jb))))
      end do
    end do
  end do
end if

do iS=1,nsym
  jS = ieor(iS-1,isym-1)+1
  if (nOrb(iS)*nOrb(jS) >= 1) then
    call DGEMM_('N','T',nOrb(iS),nOrb(jS),nOrb(iS),One,Dtmp(1+ipCM(iS)-1),nOrb(iS),rK(ipMat(jS,iS)),nOrb(jS),Zero,D(ipMat(iS,jS)), &
                nOrb(iS))
    call DGEMM_('T','N',nOrb(iS),nOrb(jS),nOrb(jS),-One,rK(ipMat(jS,iS)),nOrb(jS),Dtmp(1+ipCM(jS)-1),nOrb(jS),One,D(ipMat(iS,jS)), &
                nOrb(iS))
  end if
end do

end subroutine OITD
