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

real*8 function E2_td(FockI,rMo,loper,idisp)

use Arrays, only: G1t, G2sq
use MCLR_Data, only: nCMO, nNA, ipCM, nA
use input_mclr, only: nSym, nAsh, nIsh, nBas, ntPert
use Constants, only: Zero, Half

implicit none
integer lOper, iDisp
real*8 FockI(nCMO), rMO(*)
logical Go
integer i, j, ij, ij2, k, l, ijkl, kl2, ijkl2
integer iS, iA, iAA, iAB, jS, jA, JAA, JAB, ipF
real*8 E22
! Statement function
integer itri
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
!                                                                      *
!***********************************************************************
!                                                                      *

E22 = Zero
if (loper == 0) then
  Go = (idisp < 0)
  if (.not. Go) Go = (iand(ntpert(idisp),2**2) == 4)
  if (Go) then
    do i=1,nna
      do j=1,nna
        ij2 = i+(j-1)*nna
        ij = itri(i,j)
        do k=1,nna
          do l=1,nna
            ijkl = itri(ij,itri(k,l))
            kl2 = k+(l-1)*nna
            ijkl2 = ij2+(kl2-1)*nna**2
            E22 = E22+Half*G2sq(ijkl2)*rmo(ijkl)
          end do
        end do
      end do
    end do
  end if
  do is=1,nSym
    do iA=1,Nash(is)
      iAA = nA(is)+ia
      iAB = ia+nish(is)
      js = is
      do jA=1,nAsh(js)
        jAA = ja+nA(js)
        jAB = jA+nIsh(js)
        ipF = (iab-1)*nbas(is)+jab+ipCM(is)-1
        E22 = E22+Focki(ipf)*G1t(itri(iaa,jaa))
      end do
    end do
  end do
end if

e2_td = e22
!                                                                      *
!***********************************************************************
!                                                                      *

end function E2_td
