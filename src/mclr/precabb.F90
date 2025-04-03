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

subroutine Precabb(ib,is,js,nd,nba,rout,Temp1,ntemp,Scr,Temp2,fockti,focki,sign)
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
use MCLR_Data, only: G1t, G2t
use MCLR_Data, only: nA, nB
use input_mclr, only: nSym, nAsh, nIsh, nBas
use Constants, only: Zero, Two, Four

implicit none
integer ib, is, js, nd, nba
real*8 rout(*)
integer nTemp
real*8 Temp1(nBa,nBa), Scr(nTemp)
real*8 Temp2(nBa,nBa)
real*8 Fockti
real*8 Focki(nBa,nBa)
real*8 Sign
integer iib, jVert, ip, kS, kBB, kkB, kkC, lB, jB, ii, ij, kCC
real*8 rf, rDens, Rho

!                                                                      *
!***********************************************************************
!                                                                      *
iib = ib+nA(is)
jVert = nBas(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

ip = nTri_Elem(nd)-nTri_Elem(jVert)+1
rF = sign*Fockti
Temp2(:,:) = Zero

do kS=1,nSym
  if (nBas(js)*nash(ks) > 0) then
    do kBB=nish(ks)+1,nB(kS)
      do kCC=nish(ks)+1,kBB
        call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)
        if ((kBB > nish(ks)) .and. (kCC > nish(ks))) then
          kkB = kBB+nA(ks)-nish(ks)
          kkC = kCC+nA(ks)-Nish(ks)
          rDens = sign*Two*G2t(iTri(nTri_Elem(iib),iTri(kkb,kkc)))

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
          rDens = sign*Four*G2t(iTri(iTri(iib,kkc),iTri(kkb,iib)))
          Temp2(:,:) = Temp2(:,:)+rDens*Temp1(:,:)
        end if
      end do
    end do
  end if
end do

rho = sign*Two*G1t(nTri_Elem(iib))
do iI=nAsh(js)+nIsh(js)+1,nBas(js)
  rOut(ip) = rout(ip)-Two*rF+Rho*FockI(iI,ii)+Temp2(ii,ii)
  ip = ip+1
  do iJ=iI+1,Nbas(js)
    rOut(ip) = Rho*FockI(iI,iJ)+Temp2(ii,ij)
    ip = ip+1
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precabb
