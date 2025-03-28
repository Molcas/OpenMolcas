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

subroutine Precaaa(iC,is,js,nd,ir,rOut,nbaj,focki,fock,sign,Scr,nScr,ActInt)
!***********************************************************************
!                                                                      *
!                                        [2]                           *
!   Calculates the diagonal submatrix of E    that couple              *
!                                                                      *
!   kappa                with   kappa                for a             *
!        kactive,kactive            kactive,kactive                    *
!                                                                      *
!   single active index.                                               *
!   Used for preconditioner.                                           *
!                                                                      *
!   See Olsen,Yeager, Joergensen:                                      *
!    "Optimization and characterization of an MCSCF state"             *
!                                                                      *
!     Eq. C.12e                                                        *
!                                                                      *
!   Called by prec                                                     *
!                                                                      *
!   ib,is       :       active index for the submatrix                 *
!   js          :       symmetry of inactive,inactive                  *
!   rOut        :       Submatrix                                      *
!                                                                      *
!***********************************************************************

use Arrays, only: G1t, G2t
use MCLR_Data, only: nA
use input_mclr, only: ntAsh, nSym, nAsh, nIsh, nRS1, nRS2, nRS3
use Constants, only: Zero, One, Two, Four

implicit none
integer iC, iS, jS, nD, iR
real*8 rout(nd*(nd+1)/2)
integer nbaj
real*8 Fock(nbaj,nbaj), Focki(nbaj,nbaj)
real*8 Sign
integer nScr
real*8 Scr(nScr)
real*8 ActInt(ntAsh,ntAsh,ntAsh,ntAsh)
integer nTri, iCC, iA, iAA, jB, jBB, jjB, jD, jDD, jjD, kS, jE, jEE, jF, jFF
real*8 aecf, bedf, becf, aedf, rdbedf, rdaecf, rdaedf, rdbecf, acef, bdef, bcef, adef, rdbdef, rdacef, rdadef, rdbcef, rdbd, rdad, &
       rdac, rdbc
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)
iCC = iC+nA(iS)
!iiC = iC+nIsh(iS)
iA = iC
iAA = iCC

i = itri(iAA,iCC)
! Construct for all active orbitals first
call DCopy_(ntAsh*(ntAsh+1)/2,[Zero],0,Scr,1)
Scr(i) = One
!                                                                      *
!***********************************************************************
!                                                                      *
do jB=1,nAsh(jS) !! index B
  jBB = jB+nA(jS)
  jjB = jB+nIsh(jS)
  do jD=1,jB    !! index D
    jDD = jD+nA(jS)
    jjD = jD+nIsh(jS)
    !i = itri1(jjB,jjD)
    i = itri(jBB,jDD)
    do kS=1,nSym
      !call Coul(kS,kS,jS,jS,jB,jD,A_J,Scr)
      !call Exch(kS,jS,kS,jS,jB,jD,A_K,Scr)
      do jE=1,nAsh(ks) !! index E
        jEE = jE+nA(ks)
        do jF=1,nAsh(kS) !! index F
          jFF = jF+nA(ks)

          ! first term
          aecf = ActInt(iAA,jEE,iCC,jFF)
          bedf = ActInt(jBB,jEE,jDD,jFF)
          becf = ActInt(jBB,jEE,iCC,jFF)
          aedf = ActInt(iAA,jEE,jDD,jFF)
          rDbedf = G2t(itri(itri(jBB,jEE),itri(jDD,jFF)))
          rDaecf = G2t(itri(itri(iAA,jEE),itri(iCC,jFF)))
          rDaedf = G2t(itri(itri(iAA,jEE),itri(jDD,jFF)))
          rDbecf = G2t(itri(itri(jBB,jEE),itri(iCC,jFF)))
          Scr(i) = Scr(i)+Four*(aecf*rDbedf+bedf*rDaecf-becf*rDaedf-aedf*rDbecf)*sign

          ! second term
          acef = ActInt(iAA,iCC,jEE,jFF)
          bdef = ActInt(jBB,jDD,jEE,jFF)
          bcef = ActInt(jBB,iCC,jEE,jFF)
          adef = ActInt(iAA,jDD,jEE,jFF)
          rDbdef = G2t(itri(itri(jBB,jDD),itri(jEE,jFF)))
          rDacef = G2t(itri(itri(iAA,iCC),itri(jEE,jFF)))
          rDadef = G2t(itri(itri(iAA,jDD),itri(jEE,jFF)))
          rDbcef = G2t(itri(itri(jBB,iCC),itri(jEE,jFF)))
          Scr(i) = Scr(i)+Two*(acef*rDbdef+bdef*rDacef-bcef*rDadef-adef*rDbcef)*sign
        end do
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Work(ipG1) -> \rho_{ab}^{(1)}
i = 0 ! dummy initialize

! the remaining third and fourth terms
do jB=1,nAsh(jS)
  jBB = jB+nA(jS)
  jjB = jB+nIsh(jS)
  do jD=1,jB
    jDD = jD+nA(jS)
    jjD = jD+nIsh(jS)

    !i = itri1(jjB,jjD)
    i = itri(jBB,jDD)

    rDbd = G1t(itri(jBB,jDD))
    rDac = G1t(itri(iAA,iCC))
    rDad = G1t(itri(iAA,jDD))
    rDbc = G1t(itri(jBB,iCC))

    ! third term
    Scr(i) = Scr(i)+sign*Two*(rDbd*Focki(iA+nIsh(iS),iC+nIsh(iS))+rDac*Focki(jB+nIsh(jS),jD+nIsh(jS))- &
                              rDad*Focki(jB+nIsh(jS),iC+nIsh(iS))-rDbc*Focki(iA+nIsh(iS),jD+nIsh(jS)))
    ! fourth term
    if (iA == jD) Scr(i) = Scr(i)+sign*Two*Fock(iC+nIsh(iS),jB+nIsh(jS))
    if (jB == iC) Scr(i) = Scr(i)+sign*Two*Fock(jD+nIsh(jS),iA+nIsh(iS))
    if (jB == jD) Scr(i) = Scr(i)-sign*Two*Fock(iC+nIsh(iS),iA+nIsh(iS))
    if (iA == iC) Scr(i) = Scr(i)-sign*Two*Fock(jD+nIsh(jS),jB+nIsh(jS))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Then, Scr -> rOut
! Scr is nAsh*(nAsh+1)/2 dimension
! rOut is iOrb-nRASx space, where x is the space iC belongs to.
!nseq = 0
!do i=1,5
!  do j=1,5
!    nseq = nseq + 1
!    nseq = itri(i,j)
!    a_j(i+5*(j-1)) = scr(nseq)
!    a_j(j+5*(i-1)) = scr(nseq)
!  end do
!end do
!call sqprt(a_j,5)
!write(u6,*) 'ir = ',ir
if (iR == 1) then
  do jB=nRs1(jS)+1,nAsh(jS)
    jBB = jB+nA(jS)
    jjB = jB-nRs1(jS)+nIsh(jS)
    do jD=nRs1(jS)+1,jB
      jDD = jD+nA(jS)
      jjD = jD-nRs1(jS)+nIsh(jS)
      i = itri(jBB,jDD)
      j = itri1(jjB,jjD)
      rOut(j) = rOut(j)+Scr(i)
    end do
  end do
else if (iR == 2) then
  do jB=1,nRs1(jS)+nRs3(jS)
    jBB = jB+nA(jS)
    if (jB > nRs1(jS)) jBB = jBB+nRs2(jS)
    jjB = jB+nIsh(jS)
    do jD=1,jB
      jDD = jD+nA(jS)
      if (jD > nRs1(jS)) jDD = jDD+nRs2(jS)
      jjD = jD+nIsh(jS)
      i = itri(jBB,jDD)
      j = itri1(jjB,jjD)
      rOut(j) = rOut(j)+Scr(i)
    end do
  end do
else if (iR == 3) then
  do jB=1,nRs1(jS)+nRs2(jS)
    jBB = jB+nA(jS)
    jjB = jB+nIsh(jS)
    do jD=1,jB
      jDD = jD+nA(jS)
      jjD = jD+nIsh(jS)
      i = itri(jBB,jDD)
      j = itri1(jjB,jjD)
      rOut(j) = rOut(j)+Scr(i)
    end do
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precaaa
