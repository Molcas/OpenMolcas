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

subroutine Precabb(ib,is,js,nd,nba,rout,Temp1,ntemp,Scr,Temp2,fockti,focki,Sgn)
!***********************************************************************
!                                        [2]
!   Calculates the diagonal submatrix of E    that couple
!
!   kappa           with   kappa                for a
!        kactive,virtual        kactive,virtual
!
!   single active index.
!   Used for preconditioner.
!
!   See Olsen,Yeager, Joergensen:
!    "Optimization and characterization of an MCSCF state"
!
!   Called by prec
!
!   ib,is       :       active index for the submatrix
!   js          :       symmetry of virtual,virtual
!   rOut        :       Submatrix
!
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: G1t, G2t, nA, nB
use input_mclr, only: nAsh, nBas, nIsh, nSym
use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ib, is, js, nd, nba, nTemp
real(kind=wp), intent(inout) :: rout(*)
real(kind=wp), intent(out) :: Temp1(nBa,nBa), Scr(nTemp), Temp2(nBa,nBa)
real(kind=wp), intent(in) :: Fockti, Focki(nBa,nBa), Sgn
integer(kind=iwp) :: ii, iib, ip, jB, jVert, kBB, kCC, kkB, kkC, kS, lB
real(kind=wp) :: rf, rDens, Rho

!                                                                      *
!***********************************************************************
!                                                                      *
iib = ib+nA(is)
jVert = nBas(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

rF = Sgn*Fockti
Temp2(:,:) = Zero

do kS=1,nSym
  if (nBas(js)*nash(ks) > 0) then
    do kBB=nish(ks)+1,nB(kS)
      do kCC=nish(ks)+1,kBB
        call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)
        if ((kBB > nish(ks)) .and. (kCC > nish(ks))) then
          kkB = kBB+nA(ks)-nish(ks)
          kkC = kCC+nA(ks)-Nish(ks)
          rDens = Sgn*Two*G2t(iTri(nTri_Elem(iib),iTri(kkb,kkc)))

          if (kbb /= kcc) rDens = rDens*Two
          Temp2(:,:) = Temp2(:,:)+rDens*Temp1(:,:)
        end if
      end do
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do kS=1,nsym
  if (nBas(jS)*nAsh(kS) /= 0) then

    do lB=nish(kS)+1,nB(kS)
      kkC = nA(kS)+lB-nIsh(kS)
      do jB=nIsh(kS)+1,nB(kS)
        kkb = nA(kS)+jB-nIsh(kS)
        call EXCH(js,ks,js,ks,jb,lb,Temp1,Scr)
        if ((lB > nIsh(kS)) .and. (jB > nIsh(kS))) then
          rDens = Sgn*Four*G2t(iTri(iTri(iib,kkc),iTri(kkb,iib)))
          Temp2(:,:) = Temp2(:,:)+rDens*Temp1(:,:)
        end if
      end do
    end do
  end if
end do

ip = nTri_Elem(nd)-nTri_Elem(jVert)+1
rho = Sgn*Two*G1t(nTri_Elem(iib))
do iI=nAsh(js)+nIsh(js)+1,nBas(js)
  rOut(ip) = rout(ip)-Two*rF+Rho*FockI(iI,ii)+Temp2(ii,ii)
  rOut(ip+1:ip+nBas(jS)-iI) = Rho*FockI(iI,iI+1:nBas(jS))+Temp2(iI,iI+1:nBas(jS))
  ip = ip+nBas(jS)-iI+1
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precabb
