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

function E2(FockI,rMo,loper,idisp)

use Index_Functions, only: iTri
use MCLR_Data, only: G1t, G2t, ipCM, nA, nCMO, nNA
use input_mclr, only: nAsh, nIsh, nOrb, nSym, ntPert
use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: E2
real(kind=wp), intent(in) :: FockI(nCMO), rMO(*)
integer(kind=iwp), intent(in) :: lOper, iDisp
integer(kind=iwp) :: i, iA, iAA, iAB, ij, ijkl, ipF, iS, j, jA, jAA, jAB, jS, k, l
logical(kind=iwp) :: Go

!                                                                      *
!***********************************************************************
!                                                                      *
E2 = Zero
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
            E2 = E2+Half*G2t(ijkl)*rmo(ijkl)
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
        E2 = E2+Focki(ipf)*G1t(iTri(iaa,jaa))
      end do
    end do
  end do
end if

!                                                                      *
!***********************************************************************
!                                                                      *

end function E2
