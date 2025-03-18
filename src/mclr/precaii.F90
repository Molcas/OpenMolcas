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

subroutine Precaii(iB,is,js,nd,ir,rOut,nbai,nbaj,fockii,fockai,fockti,focki,focka,fock,sign,A_J,A_K,Scr,nScr)
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

use Arrays, only: G1t, G2t
use MCLR_Data, only: nA
use input_mclr, only: nSym, nAsh, nIsh, nBas

implicit none
integer iB, is, js, nd, ir
real*8 rout(nd*(nd+1)/2)
integer nbai, nbaj
real*8 fockii, fockai, fockti
real*8 Fock(nbaj,nbaj), FockA(nBaj,nBaj), Focki(nbaj,nbaj)
integer nScr
real*8 A_J(nScr), A_K(nScr), Scr(nScr)
real*8 sign
integer nTri, iBB, iiB, jA, jB, kS, jC, jCC, jjC, jD, jDD, jjD, iCD, iC, iCC, iiC, iCB, iBC
real*8 rDens1, rDens2, CDij, CiDj, rDens, CiBj, BiCj, BCij, rFock
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)
iBB = ib+nA(is)
iiB = ib+nish(is)
!                                                                      *
!***********************************************************************
!                                                                      *
do jA=1,nIsh(jS)
  do jB=1,jA

    i = itri1(ja,jb)

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

          rDens1 = sign*G2t(itri(itri(jCC,jDD),itri(iBB,iBB)))

          ! gamma(bdcb)

          rDens2 = sign*G2t(itri(itri(iBB,jDD),itri(jCC,iBB)))

          ! (cd|ij)

          icd = (jjD-1)*nBas(kS)+jjC
          cdij = A_J(icd)
          cidj = A_K(icd)
          ! is the coefficient opposite?
          rout(i) = rout(i)+2.0d0*(rDens2*cidj+2.0d0*rDens1*cdij)
          !rout(i) = rout(i)+2.0d0*(rDens1*cdij+2.0d0*rDens2*cidj)
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

  rDens = sign*(-G1t(itri(iCC,iBB)))
  if (iCC == iBB) rdens = rdens+sign
  rDens = 2.0d0*rDens

  call Coul(jS,jS,iS,iS,iiB,iiC,A_J,Scr)
  call Exch(jS,iS,jS,iS,iiB,iiC,A_K,Scr)

  do jA=1,nIsh(jS)
    do jB=1,jA
      i = itri1(jA,jB)

      ! (ci|bj)
      ! (bi|cj)

      icb = (jB-1)*nBas(jS)+jA
      ibc = (jA-1)*nBas(jS)+jB
      cibj = A_K(icb)
      bicj = A_K(ibc)

      ! (bc|ij)

      bcij = A_J(ibc)

      rout(i) = rout(i)+rdens*(7.0d0*cibj-sign*bicj-sign*2.0d0*bcij)

    end do
  end do

end do
!                                                                      *
!***********************************************************************
!                                                                      *
rFock = sign*2.0d0*Fockii+sign*2.0d0*Fockai-sign*Fockti
rdens = sign*2.0d0*G1t(itri(ibb,ibb))
i = 0 ! dummy initialize

do jA=1,nIsh(jS)
  do jB=1,jA

    i = itri1(ja,jb)

    rout(i) = rout(i)-sign*4.0d0*(Focka(jA,jB)+Focki(jA,jB))+rdens*Focki(ja,jb)
  end do
  rout(i) = rout(i)+2.0d0*rfock
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(ir)
  call Unused_integer(nbai)
  call Unused_real_array(fock)
end if

end subroutine Precaii
