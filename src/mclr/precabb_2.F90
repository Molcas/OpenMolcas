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

subroutine Precabb_2(ib,is,js,nd,no,rout,Temp1,ntemp,Scr,Temp2,fockti,focki,Sgn)
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
use input_mclr, only: nAsh, nIsh, nOrb, nSym
use Constants, only: Zero, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ib, is, js, nd, no, nTemp
real(kind=wp), intent(inout) :: rout(*)
real(kind=wp), intent(out) :: Temp1(no,no), Scr(nTemp), Temp2(nO,nO)
real(kind=wp), intent(in) :: Fockti, Focki(no,no), Sgn
integer(kind=iwp) :: ii, iib, ijkl, ip, jB, jVert, kBB, kCC, kkB, kkC, kS, lB
real(kind=wp) :: rDens1, rDens2, rf, Rho

!                                                                      *
!***********************************************************************
!                                                                      *
iib = ib+nA(is)
jVert = nOrb(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

rF = Sgn*Fockti
Temp2(:,:) = Zero

do kS=1,nSym
  if (nOrb(js)*nash(ks) > 0) then

    do kBB=nish(ks)+1,nB(kS)
      do kCC=nish(ks)+1,kBB
        call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)

        if ((kBB > nish(ks)) .and. (kCC > nish(ks))) then
          kkB = kBB+nA(ks)-nish(ks)
          kkC = kCC+nA(ks)-Nish(ks)
          rDens1 = Sgn*Two*G2t(iTri(nTri_Elem(iib),iTri(kkb,kkc)))

          if (kbb /= kcc) rdens1 = rdens1*Two

          Temp2(:,:) = Temp2(:,:)+rdens1*Temp1(:,:)

        end if
      end do
    end do
  end if
end do

do Ks=1,nsym
  ijkl = nOrb(js)*nash(ks)
  if (ijkl /= 0) then

    do LB=nish(ks)+1,nB(KS)
      kkc = nA(ks)+lb-nish(ks)
      do JB=nish(ks)+1,nB(KS)
        kkb = nA(ks)+jb-nish(ks)
        call EXCH(js,ks,js,ks,jb,lb,Temp1,Scr)
        if ((LB > nISH(ks)) .and. (jb > nish(ks))) then
          rDens2 = Sgn*Four*G2t(iTri(iTri(iib,kkc),iTri(kkb,iib)))
          Temp2(:,:) = Temp2(:,:)+rDens2*Temp1(:,:)
        end if
      end do
    end do
  end if
end do

ip = nTri_Elem(nd)-nTri_Elem(jVert)+1
rho = Sgn*Two*G1t(nTri_Elem(iib))
do iI=nAsh(js)+nIsh(js)+1,nOrb(js)
  rOut(ip) = rout(ip)-Two*rF+Rho*FockI(iI,ii)+Temp2(ii,ii)
  rOut(ip+1:ip+nOrb(jS)-iI) = Rho*FockI(iI,iI+1:nOrb(jS))+Temp2(iI,iI+1:nOrb(jS))
  ip = ip+nOrb(jS)-iI+1
end do

end subroutine Precabb_2
