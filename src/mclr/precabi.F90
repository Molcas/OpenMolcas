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

subroutine Precabi(ib,is,js,ir,nd,rOut,nba,focki,focka,fock,sign,A_J,A_K,Scr,nScr)
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

use Arrays, only: G1t, G2t
use MCLR_Data, only: nA
use input_mclr, only: nSym, nAsh, nIsh, nOrb, nBas

implicit none
integer ib, is, js, ir, nd
real*8 rOut(*)
integer nba
real*8 Fock(nba,nba), Focki(nba,nba), FockA(nba,nba)
real*8 Sign
integer nScr
real*8 A_J(nScr), A_K(nScr), Scr(nScr)
integer nTri, jVert, nO, iAA, itAA, kS, kA, kAA, kkA, lA, lAA, llA, jB, ip, iVJ
real*8 Fact, Fact1, Fact2
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)

jVert = nOrb(js)-nIsh(js)-nAsh(js)
if (jVert == 0) return

nO = nAsh(js)+nIsh(js)
iAA = nA(is)+ib
itAA = itri(iAA,iAA)
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
        ip = itri1(jB,nd-jVert+1)

        Fact1 = -2.0d0*G2t(itri(itAA,itri(kAA,lAA)))
        Fact2 = -4.0d0*G2t(itri(itri(iAA,kAA),itri(iAA,lAA)))

        if (kaa == iaa) Fact2 = Fact2+8.0d0*G1t(itri(iAA,lAA))
        if (laa == iaa) Fact1 = Fact1-2.0d0*G1t(itri(iAA,kAA))
        if (laa == iaa) Fact2 = Fact2-2.0d0*G1t(itri(iAA,kAA))

        ivj = (jB-1)*nBas(jS)+no+1
        call DaXpY_(jVert,Sign*Fact1,A_J(ivj),1,rout(ip),1) ! ????
        call DaXpY_(jVert,Sign*Fact2,A_K(ivj),1,rout(ip),1)

      end do

    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do jB=1,nIsh(jS)
  ip = itri1(jB,nd-jVert+1)
  Fact = (2.0d0-2.0d0*G1t(itAA))
  call DaxPy_(jVert,Sign*Fact,FockI(nO+1,jB),1,rOut(ip),1)
  Fact = 2.0d0
  call DaxPy_(jVert,Sign*Fact,FockA(nO+1,jB),1,rOut(ip),1)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(ir)
  call Unused_real_array(fock)
end if

end subroutine Precabi
