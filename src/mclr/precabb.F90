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

subroutine Precabb(ib,is,js,nd,nba,rout,Temp1,ntemp,Scr,Temp2,fockti,focki,focka,fock,sign)
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
use input_mclr, only: nSym, nAsh, nIsh, nBas
use Constants, only: Zero, Two, Four

implicit none
integer ib, is, js, nd, nba
real*8 rout(*)
integer nTemp
real*8 Temp1(nTemp), Scr(nTemp)
real*8 Temp2(nBa,nBa)
real*8 Fockti
real*8 Focki(nBa,nBa), Focka(nBa,nBa), Fock(nBa,nBa)
real*8 Sign
integer nTri, iib, jVert, i2, ip, kS, kBB, kkB, kkC, lB, jB, ii, ij, kCC
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
jVert = nBas(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

i2 = nD-jVert+1
ip = iTri1(i2,i2)
rF = sign*Fockti
call dcopy_(nBa**2,[Zero],0,Temp2,1)

do kS=1,nSym
  if (nBas(js)*nash(ks) > 0) then
    do kBB=nish(ks)+1,nB(kS)
      do kCC=nish(ks)+1,kBB
        call COUL(jS,jS,kS,kS,kbb,kcc,Temp1,Scr)
        if ((kBB > nish(ks)) .and. (kCC > nish(ks))) then
          kkB = kBB+nA(ks)-nish(ks)
          kkC = kCC+nA(ks)-Nish(ks)
          rDens1 = sign*Two*G2t(itri(itri(iib,iib),itri(kkb,kkc)))

          if (kbb /= kcc) rdens1 = rdens1*Two
          call DaxPy_(nBa**2,rdens1,Temp1,1,Temp2,1)
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
          rDens2 = sign*Four*G2t(itri(itri(iib,kkc),itri(kkb,iib)))
          call DaXpY_(nBa**2,rDens2,Temp1,1,Temp2,1)
        end if
      end do
    end do
  end if
end do

rho = sign*Two*G1t(itri(iib,iib))
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
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(focka)
  call Unused_real_array(fock)
end if

end subroutine Precabb
