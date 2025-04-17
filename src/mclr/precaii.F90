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

subroutine Precaii(iB,is,js,nd,rOut,nbaj,fockii,fockai,fockti,focki,focka,Sgn,A_J,A_K,Scr,nScr)
!***********************************************************************
!                                                                      *
!                                        [2]                           *
!   Calculates the diagonal submatrix of E    that couple              *
!                                                                      *
!   kappa                with   kappa                for a             *
!        kactive,inactive            kactive,inactive                  *
!                                                                      *
!   single active index.                                               *
!   Used for preconditioner.                                           *
!                                                                      *
!   See Olsen,Yeager, Joergensen:                                      *
!    "Optimization and characterization of an MCSCF state"             *
!                                                                      *
!     Eq. C.12a                                                        *
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
use input_mclr, only: nAsh, nBas, nIsh, nSym
use Constants, only: Two, Four, Seven
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iB, is, js, nd, nbaj, nScr
real(kind=wp), intent(inout) :: rout(nTri_Elem(nd))
real(kind=wp), intent(in) :: fockii, fockai, fockti, Focki(nbaj,nbaj), FockA(nBaj,nBaj), Sgn
real(kind=wp), intent(out) :: A_J(nScr), A_K(nScr), Scr(nScr)
integer(kind=iwp) :: i, iBB, iBC, iC, iCB, iCC, iCD, iiB, iiC, jA, jB, jC, jCC, jD, jDD, jjC, jjD, kS, nTri
real(kind=wp) :: BCij, BiCj, CDij, CiBj, CiDj, rDens, rDens1, rDens2, rFock

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = nTri_Elem(nd)
iBB = ib+nA(is)
iiB = ib+nish(is)
!                                                                      *
!***********************************************************************
!                                                                      *
do jA=1,nIsh(jS)
  do jB=1,jA

    i = nTri-iTri(nd-ja+1,nd-jb+1)+1

    do kS=1,nSym

      call Coul(kS,kS,jS,jS,jA,jB,A_J,Scr)
      call Exch(kS,jS,kS,jS,jA,jB,A_K,Scr)

      do jC=1,nAsh(ks)
        jCC = jC+nA(ks)
        jjC = jC+nIsh(ks)
        do jD=1,nAsh(kS)
          jDD = jD+nA(ks)
          jjD = jD+nIsh(ks)

          ! gamma(cdbb)=gamma(bbcd)

          rDens1 = Sgn*G2t(iTri(iTri(jCC,jDD),nTri_Elem(iBB)))

          ! gamma(bdcb)

          rDens2 = Sgn*G2t(iTri(iTri(iBB,jDD),iTri(jCC,iBB)))

          ! (cd|ij)

          icd = (jjD-1)*nBas(kS)+jjC
          cdij = A_J(icd)
          cidj = A_K(icd)
          ! is the coefficient opposite?
          rout(i) = rout(i)+Two*(rDens2*cidj+Two*rDens1*cdij)
          !rout(i) = rout(i)+Two*(rDens1*cdij+Two*rDens2*cidj)
        end do
      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do iC=1,nAsh(is)
  iCC = ic+nA(is)
  iiC = iC+nIsh(is)

  ! 2*(delta(bc)-D(bc))

  rDens = Sgn*(-G1t(iTri(iCC,iBB)))
  if (iCC == iBB) rdens = rdens+Sgn
  rDens = Two*rDens

  call Coul(jS,jS,iS,iS,iiB,iiC,A_J,Scr)
  call Exch(jS,iS,jS,iS,iiB,iiC,A_K,Scr)

  do jA=1,nIsh(jS)
    do jB=1,jA
      i = nTri-iTri(nd-jA+1,nd-jB+1)+1

      ! (ci|bj)
      ! (bi|cj)

      icb = (jB-1)*nBas(jS)+jA
      ibc = (jA-1)*nBas(jS)+jB
      cibj = A_K(icb)
      bicj = A_K(ibc)

      ! (bc|ij)

      bcij = A_J(ibc)

      rout(i) = rout(i)+rdens*(Seven*cibj-Sgn*(bicj+Two*bcij))

    end do
  end do

end do
!                                                                      *
!***********************************************************************
!                                                                      *
rFock = Sgn*(Two*(Fockii+Fockai)-Fockti)
rdens = Sgn*Two*G1t(nTri_Elem(ibb))
i = 0 ! dummy initialize

do jA=1,nIsh(jS)
  do jB=1,jA

    i = nTri-iTri(nd-ja+1,nd-jb+1)+1

    rout(i) = rout(i)-Sgn*Four*(Focka(jA,jB)+Focki(jA,jB))+rdens*Focki(ja,jb)
  end do
  rout(i) = rout(i)+Two*rfock
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precaii
