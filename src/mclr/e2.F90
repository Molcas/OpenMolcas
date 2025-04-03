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

real*8 function E2(FockI,rMo,loper,idisp)

use Index_Functions, only: iTri
use MCLR_Data, only: G1t, G2t
use MCLR_Data, only: nCMO, nNA, ipCM, nA
use input_mclr, only: nSym, nAsh, nIsh, nOrb, ntPert
use Constants, only: Zero, Half

implicit none
integer lOper, iDisp
real*8 FockI(nCMO), rMO(*)
logical Go
real*8 E22
integer i, j, ij, k, l, ijkl, iS, jS, iA, jA, iAA, iAB, jAA, jAB, ipF

!                                                                      *
!***********************************************************************
!                                                                      *
E22 = Zero
if (loper == 0) then
  Go = (iDisp < 0)
  if (.not. Go) Go = btest(ntpert(idisp),2)
  if (Go) then
    do i=1,nna
      do j=1,nna
        ij = iTri(i,j)
        do k=1,nna
          do l=1,nna
            ijkl = iTri(ij,iTri(k,l))
            E22 = E22+Half*G2t(ijkl)*rmo(ijkl)
          end do
        end do
      end do
    end do
  end if
  do is=1,nSym
    do iA=1,nAsh(is)
      iAA = nA(iS)+ia
      iAB = ia+nIsh(iS)
      js = is
      do jA=1,nAsh(js)
        jAA = ja+nA(js)
        jAB = jA+nIsh(js)
        ipF = (iab-1)*norb(is)+jab+ipCM(is)-1
        E22 = E22+Focki(ipf)*G1t(iTri(iaa,jaa))
      end do
    end do
  end do
end if

e2 = e22
!                                                                      *
!***********************************************************************
!                                                                      *

end function E2
