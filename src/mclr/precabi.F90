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
! Copyright (C) 1996, Anders Bernhardsson                              *
!***********************************************************************

subroutine Precabi(ib,is,js,nd,rOut,nba,focki,focka,Sgn,A_J,A_K,Scr,nScr)
!***********************************************************************
!                                          [2]                         *
!     Calculates the diagonal submatrix of E    that couple            *
!                                                                      *
!     kappa           with   kappa                for a                *
!          kactive,virtual        kactive,inactive                     *
!                                                                      *
!     single active index.                                             *
!     Used for preconditioner.                                         *
!                                                                      *
!     See Olsen,Yeager, Joergensen:                                    *
!      "Optimization and characterization of an MCSCF state"           *
!                                                                      *
!     Called by prec                                                   *
!                                                                      *
!     ib,is         :        active index for the submatrix            *
!     js            :        symmetry of inact,virtual                 *
!     rOut          :         Submatrix                                *
!                                                                      *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: G1t, G2t, nA
use input_mclr, only: nAsh, nBas, nIsh, nOrb, nSym
use Constants, only: Two, Four, Eight
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ib, is, js, nd, nba, nScr
real(kind=wp), intent(inout) :: rOut(*)
real(kind=wp), intent(in) :: Focki(nba,nba), FockA(nba,nba), Sgn
real(kind=wp), intent(out) :: A_J(nScr), A_K(nScr), Scr(nScr)
integer(kind=iwp) :: iAA, ip, itAA, iVJ, jB, jVert, kA, kAA, kkA, kS, lA, lAA, llA, nO, nTri
real(kind=wp) :: Fact, Fact1, Fact2

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = nTri_Elem(nd)

jVert = nOrb(js)-nIsh(js)-nAsh(js)
if (jVert == 0) return

nO = nAsh(js)+nIsh(js)
iAA = nA(is)+ib
itAA = nTri_Elem(iAA)
!                                                                      *
!***********************************************************************
!                                                                      *
do kS=1,nSym
  do kA=1,nAsh(kS)
    kAA = kA+nA(ks)
    kkA = kA+nIsh(ks)
    do lA=1,nAsh(kS)
      lAA = lA+nA(ks)
      llA = lA+nIsh(ks)

      call Coul(jS,jS,kS,kS,kkA,llA,A_J,Scr)
      call Exch(jS,kS,jS,kS,kkA,llA,A_K,Scr)

      do jB=1,nIsh(jS)
        ip = nTri-iTri(nd-jB+1,jVert)

        Fact1 = -Two*G2t(iTri(itAA,iTri(kAA,lAA)))
        Fact2 = -Four*G2t(iTri(iTri(iAA,kAA),iTri(iAA,lAA)))

        if (kaa == iaa) Fact2 = Fact2+Eight*G1t(iTri(iAA,lAA))
        if (laa == iaa) Fact1 = Fact1-Two*G1t(iTri(iAA,kAA))
        if (laa == iaa) Fact2 = Fact2-Two*G1t(iTri(iAA,kAA))

        ivj = (jB-1)*nBas(jS)+no
        rout(ip+1:ip+jVert) = rout(ip+1:ip+jVert)+Sgn*(Fact1*A_J(ivj+1:ivj+jVert)+Fact2*A_K(ivj+1:ivj+jVert))

      end do

    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do jB=1,nIsh(jS)
  ip = nTri-iTri(nd-jB+1,jVert)
  Fact = (Two-Two*G1t(itAA))
  rOut(ip+1:ip+jVert) = rOut(ip+1:ip+jVert)+Sgn*(Fact*FockI(nO+1:nO+jVert,jB)+Two*FockA(nO+1:nO+jVert,jB))
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precabi
