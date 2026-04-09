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

subroutine Preciaa(iB,iS,jS,nd,rOut,nbaj,fockii,fockai,focki,focka,fock,Sgn,A_J,A_K,Scr,nScr)
!***********************************************************************
!     Change Fock(i) ne Fock(j)
!                                           [2]
!     Calculates the diagonal submatrix of E    that couple
!
!     kappa               with   kappa                for a
!          kinactive,active           kinactive,active
!
!     single inactive index.
!     Used for preconditioner.
!
!     See Olsen,Yeager, Joergensen:
!      "Optimization and characterization of an MCSCF state"
!
!     Called by prec
!
!     ib,is       :       inactive index for the submatrix
!     js          :       symmetry of active,active
!     rOut        :       Submatrix
!
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem
use MCLR_Data, only: G1t, G2t, nA
use input_mclr, only: nAsh, nBas, nIsh, nSym
use Constants, only: Two, Three, Four
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iB, iS, jS, nd, nbaj, nScr
real(kind=wp), intent(inout) :: rout(*)
real(kind=wp), intent(in) :: fockii, fockai, focki(nbaj,nbaj), focka(nbaj,nbaj), fock(nbaj,nbaj), Sgn
real(kind=wp), intent(out) :: A_J(nScr), A_K(nScr), Scr(nScr)
integer(kind=iwp) :: i, iAC, iBC, ip1, ip2, jA, jAA, jB, jBB, jC, jCC, jD, jDD, jjA, jjB, jjC, jjD, kS, nTri
real(kind=wp) :: AABB, ABAB, ABCB, ACBB, BBCB, BCBB, rDens, rDens1, rDens2, rFock

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = nTri_Elem(nd)
!                                                                      *
!***********************************************************************
!                                                                      *
do kS=1,nSym

  if (nAsh(kS) /= 0) call Coul(kS,kS,iS,iS,iB,iB,A_J,Scr)

  do jC=1,nAsh(kS)
    jjC = JC+nA(kS)
    jCC = jC+nIsh(kS)

    call Coul(kS,iS,kS,iS,jCC,iB,A_K,Scr)

    do jD=1,nAsh(kS)
      jjD = JD+nA(kS)
      jDD = jD+nIsh(kS)

      ip1 = (jDD-1)*nBas(kS)+jCC
      ip2 = (iB-1)*nBas(kS)+jDD
      aabb = A_J(ip1)
      abab = A_K(ip2)

      do jA=1,nAsh(jS)
        jjA = jA+nA(jS)
        do jB=1,jA
          jjB = jB+nA(jS)
          i = nTri-iTri(nd-jA+1,nd-jB+1)+1

          rDens1 = Sgn*G2t(iTri(iTri(jjC,jjD),iTri(jjB,jjA)))
          rDens2 = Sgn*G2t(iTri(iTri(jjB,jjD),iTri(jjC,jjA)))

          rout(i) = rout(i)+Two*rDens1*aabb+Four*rDens2*abab

        end do
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (nAsh(jS) /= 0) then
  call Exch(jS,iS,jS,iS,iB,iB,A_K,Scr)
  call Coul(jS,jS,iS,iS,iB,iB,A_J,Scr)
end if

do jA=1,nAsh(jS)
  jAA = jA+nA(jS)
  jjA = jA+nIsh(jS)

  do jB=1,jA
    jBB = jB+nA(jS)
    jjB = jB+nIsh(jS)
    i = nTri-iTri(nd-jA+1,nd-jB+1)+1

    do jC=1,nAsh(jS)
      jCC = jC+nA(jS)
      jjC = jC+nIsh(jS)
      iBC = (jjC-1)*nBas(jS)+jjB
      BCbb = A_J(iBC)
      BbCb = A_K(iBC)
      iAC = (jjC-1)*nBas(jS)+jjA
      ACbb = A_J(iAC)
      AbCb = A_K(iAC)

      rDens1 = -Sgn*G1t(iTri(jAA,jCC))
      rDens2 = -Sgn*G1t(iTri(jBB,jCC))
      if (jAA == jCC) rDens1 = rdens1+Sgn
      if (jBB == jCC) rDens2 = rdens2+Sgn

      rout(i) = rout(i)+Two*rdens1*(Three*BbCb-BCbb)+Two*rdens2*(Three*AbCb-ACbb)

    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
rFock = Sgn*(Fockii+FockAi)
i = 0 ! dummy initialize
do jA=1,nAsh(jS)
  jAA = jA+nA(jS)
  jjA = jA+nIsh(js)
  do jB=1,JA
    jBB = jB+nA(jS)
    jjB = jB+nIsh(js)
    i = nTri-iTri(nd-jA+1,nd-jB+1)+1
    rDens = G1t(iTri(jbb,jAA))
    rout(i) = rout(i)+Sgn*(Two*rdens*Fockii+Two*(Two*Focki(jjA,jjB)+Two*FockA(jjA,jjB)-Fock(jjB,jjA)))
  end do
  rout(i) = rout(i)-Four*rFock
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Preciaa
