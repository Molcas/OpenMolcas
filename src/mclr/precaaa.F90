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

subroutine Precaaa(iC,is,js,nd,ir,rOut,nbaj,focki,fock,Sgn,Scr,nScr,ActInt)
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

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: G1t, G2t, nA
use input_mclr, only: nAsh, nIsh, nRS1, nRS2, nRS3, nSym, ntAsh
use Constants, only: One, Two, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iC, iS, jS, nD, iR, nbaj, nScr
real(kind=wp), intent(inout) :: rout(nTri_Elem(nd))
real(kind=wp), intent(in) :: Focki(nbaj,nbaj), Fock(nbaj,nbaj), Sgn, ActInt(ntAsh,ntAsh,ntAsh,ntAsh)
real(kind=wp), intent(out) :: Scr(nScr)
integer(kind=iwp) :: i, iA, iAA, iCC, j, jB, jBB, jD, jDD, jE, jEE, jF, jFF, jjB, jjD, kS, nTri
real(kind=wp) :: acef, adef, aecf, aedf, bcef, bdef, becf, bedf, rdac, rdacef, rdad, rdadef, rdaecf, rdaedf, rdbc, rdbcef, rdbd, &
                 rdbdef, rdbecf, rdbedf

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = nTri_Elem(nd)
iCC = iC+nA(iS)
!iiC = iC+nIsh(iS)
iA = iC
iAA = iCC

i = iTri(iAA,iCC)
! Construct for all active orbitals first
Scr(1:nTri_Elem(ntAsh)) = 0
Scr(i) = One
!                                                                      *
!***********************************************************************
!                                                                      *
do jB=1,nAsh(jS) !! index B
  jBB = jB+nA(jS)
  !jjB = nd-(jB+nIsh(jS))+1
  do jD=1,jB    !! index D
    jDD = jD+nA(jS)
    !jjD = nd-(jD+nIsh(jS))+1
    !i = nTri-iTri(jjB,jjD)+1
    i = iTri(jBB,jDD)
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
          rDbedf = G2t(iTri(iTri(jBB,jEE),iTri(jDD,jFF)))
          rDaecf = G2t(iTri(iTri(iAA,jEE),iTri(iCC,jFF)))
          rDaedf = G2t(iTri(iTri(iAA,jEE),iTri(jDD,jFF)))
          rDbecf = G2t(iTri(iTri(jBB,jEE),iTri(iCC,jFF)))
          Scr(i) = Scr(i)+Four*(aecf*rDbedf+bedf*rDaecf-becf*rDaedf-aedf*rDbecf)*Sgn

          ! second term
          acef = ActInt(iAA,iCC,jEE,jFF)
          bdef = ActInt(jBB,jDD,jEE,jFF)
          bcef = ActInt(jBB,iCC,jEE,jFF)
          adef = ActInt(iAA,jDD,jEE,jFF)
          rDbdef = G2t(iTri(iTri(jBB,jDD),iTri(jEE,jFF)))
          rDacef = G2t(iTri(iTri(iAA,iCC),iTri(jEE,jFF)))
          rDadef = G2t(iTri(iTri(iAA,jDD),iTri(jEE,jFF)))
          rDbcef = G2t(iTri(iTri(jBB,iCC),iTri(jEE,jFF)))
          Scr(i) = Scr(i)+Two*(acef*rDbdef+bdef*rDacef-bcef*rDadef-adef*rDbcef)*Sgn
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
  !jjB = nd-(jB+nIsh(jS))+1
  do jD=1,jB
    jDD = jD+nA(jS)
    !jjD = nd-(jD+nIsh(jS))+1

    !i = nTri-iTri(jjB,jjD)+1
    i = iTri(jBB,jDD)

    rDbd = G1t(iTri(jBB,jDD))
    rDac = G1t(iTri(iAA,iCC))
    rDad = G1t(iTri(iAA,jDD))
    rDbc = G1t(iTri(jBB,iCC))

    ! third term
    Scr(i) = Scr(i)+Sgn*Two*(rDbd*Focki(iA+nIsh(iS),iC+nIsh(iS))+rDac*Focki(jB+nIsh(jS),jD+nIsh(jS))- &
                             rDad*Focki(jB+nIsh(jS),iC+nIsh(iS))-rDbc*Focki(iA+nIsh(iS),jD+nIsh(jS)))
    ! fourth term
    if (iA == jD) Scr(i) = Scr(i)+Sgn*Two*Fock(iC+nIsh(iS),jB+nIsh(jS))
    if (jB == iC) Scr(i) = Scr(i)+Sgn*Two*Fock(jD+nIsh(jS),iA+nIsh(iS))
    if (jB == jD) Scr(i) = Scr(i)-Sgn*Two*Fock(iC+nIsh(iS),iA+nIsh(iS))
    if (iA == iC) Scr(i) = Scr(i)-Sgn*Two*Fock(jD+nIsh(jS),jB+nIsh(jS))
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
!    nseq = nseq+1
!    nseq = iTri(i,j)
!    a_j(i+5*(j-1)) = scr(nseq)
!    a_j(j+5*(i-1)) = scr(nseq)
!  end do
!end do
!call sqprt(a_j,5)
!write(u6,*) 'ir = ',ir
select case (iR)
  case (1)
    do jB=nRs1(jS)+1,nAsh(jS)
      jBB = jB+nA(jS)
      jjB = nd-(jB-nRs1(jS)+nIsh(jS))+1
      do jD=nRs1(jS)+1,jB
        jDD = jD+nA(jS)
        jjD = nd-(jD-nRs1(jS)+nIsh(jS))+1
        i = iTri(jBB,jDD)
        j = nTri-iTri(jjB,jjD)+1
        rOut(j) = rOut(j)+Scr(i)
      end do
    end do
  case (2)
    do jB=1,nRs1(jS)+nRs3(jS)
      jBB = jB+nA(jS)
      if (jB > nRs1(jS)) jBB = jBB+nRs2(jS)
      jjB = nd-(jB+nIsh(jS))+1
      do jD=1,jB
        jDD = jD+nA(jS)
        if (jD > nRs1(jS)) jDD = jDD+nRs2(jS)
        jjD = nd-(jD+nIsh(jS))+1
        i = iTri(jBB,jDD)
        j = nTri-iTri(jjB,jjD)+1
        rOut(j) = rOut(j)+Scr(i)
      end do
    end do
  case (3)
    do jB=1,nRs1(jS)+nRs2(jS)
      jBB = jB+nA(jS)
      jjB = nd-(jB+nIsh(jS))+1
      do jD=1,jB
        jDD = jD+nA(jS)
        jjD = nd-(jD+nIsh(jS))+1
        i = iTri(jBB,jDD)
        j = nTri-iTri(jjB,jjD)+1
        rOut(j) = rOut(j)+Scr(i)
      end do
    end do
end select
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precaaa
