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

subroutine Precabb_2(ib,is,js,nd,nba,no,rout,Temp1,ntemp,Scr,Temp2,fockti,focki,sign)
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

use Arrays, only: G1t, G2t
use MCLR_Data, only: nA, nB
use input_mclr, only: nSym, nAsh, nIsh, nOrb
use Constants, only: Zero, Two, Four

implicit none
integer ib, is, js, nd, nba, no
real*8 rout(*)
integer nTemp
real*8 Temp1(nTemp), Temp2(nO,nO), Scr(nTemp)
real*8 Fockti
real*8 Focki(no,no)
real*8 Sign
integer nTri, iib, jVert, i2, ip, kS, kBB, ipT, kkB, kkC, ijkl, lB, jB, ii, ij, kCC
real*8 rf, rDens1, rDens2, Rho
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)

iib = ib+nA(is)
jVert = nOrb(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

i2 = nD-jVert+1
ip = iTri1(i2,i2)
rF = sign*Fockti
call dcopy_(nBa**2,[Zero],0,Temp2,1)

do kS=1,nSym
  if (nOrb(js)*nash(ks) > 0) then

    do kBB=nish(ks)+1,nB(kS)
      do kCC=nish(ks)+1,kBB
        call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)
        ipT = 1

        if ((kBB > nish(ks)) .and. (kCC > nish(ks))) then
          kkB = kBB+nA(ks)-nish(ks)
          kkC = kCC+nA(ks)-Nish(ks)
          rDens1 = sign*Two*G2t(itri(itri(iib,iib),itri(kkb,kkc)))

          if (kbb /= kcc) rdens1 = rdens1*Two

          call DaxPy_(nO**2,rdens1,Temp1,1,Temp2,1)

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
        ipT = 1
        if ((LB > nISH(ks)) .and. (jb > nish(ks))) then
          rDens2 = sign*Four*G2t(itri(itri(iib,kkc),itri(kkb,iib)))
          call DaXpY_(nO**2,rDens2,Temp1(ipT),1,Temp2,1)
        end if
      end do
    end do
  end if
end do

rho = sign*Two*G1t(itri(iib,iib))
do iI=nAsh(js)+nIsh(js)+1,nOrb(js)
  rOut(ip) = rout(ip)-Two*rF+Rho*FockI(iI,ii)+Temp2(ii,ii)
  ip = ip+1
  do iJ=iI+1,NOrb(js)
    rOut(ip) = Rho*FockI(iI,iJ)+Temp2(ii,ij)
    ip = ip+1
  end do
end do

end subroutine Precabb_2
