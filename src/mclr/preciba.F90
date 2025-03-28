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

subroutine Preciba(iB,iS,jS,nd,rOut,nba,focki,focka,fock,sign,A_J,A_K,Scr,nScr)
!***********************************************************************
!                                          [2]                         *
!     Calculates the diagonal submatrix of E    that couple            *
!                                                                      *
!     kappa           with   kappa                for a                *
!          kinactive,virtual        kinactive,active                   *
!                                                                      *
!     single inactive index.                                           *
!     Used for preconditioner.                                         *
!                                                                      *
!     See Olsen,Yeager, Joergensen:                                    *
!      "Optimization and characterization of an MCSCF state"           *
!                                                                      *
!     Called by prec                                                   *
!                                                                      *
!     ib,is       :       inactive index for the submatrix             *
!     js          :       symmetry of virtual,active                   *
!     rOut        :       Submatrix                                    *
!                                                                      *
!***********************************************************************

use MCLR_Data, only: G1t
use MCLR_Data, only: nA
use input_mclr, only: nAsh, nIsh, nBas, nOrb
use Constants, only: Two, Four, Six

implicit none
integer iB, iS, jS, nd
real*8 rOut(*)
integer nba
real*8 Fock(nba,nba), Focki(nba,nba), FockA(nba,nba)
real*8 Sign
integer nScr
real*8 A_J(nScr), A_K(nScr), Scr(nScr)
integer nTri, jVert, nO, jA, ip, jB, jBB, iVB
real*8 rDens
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)
jVert = nOrb(jS)-nAsh(jS)-nIsh(jS)
nO = nAsh(jS)+nIsh(jS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Get block J^(iB,iB)
call Coul(jS,jS,iS,iS,iB,iB,A_J,Scr)
! Get block K^(iB,iB)
call Exch(jS,iS,jS,iS,iB,iB,A_K,Scr)

do jA=1,nAsh(jS)
  ip = itri1(ja,nd-jVert+1)
  do jB=1,nAsh(jS)
    jBB = jB+nIsh(jS)
    ! Get D_(ja,jb)
    rDens = -sign*G1t((iTri(jA+nA(jS),jB+nA(jS))))
    if (jA == jB) rDens = rdens+sign*Two

    ivB = (jBB-1)*nBas(jS)+nO+1
    call DaXpY_(jVert,Six*rDens,A_K(ivB),1,rOut(ip),1) ! ????
    call DaXpY_(jVert,-Two*rDens,A_J(ivB),1,rOut(ip),1)
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do jA=1,nAsh(js)
  ip = iTri1(ja,nAsh(js)+1)
  call DaXpY_(jVert,sign*Four,Focki(nO+1,ja+nIsh(js)),1,rout(ip),1)
  call DaXpY_(jVert,sign*Four,FockA(nO+1,ja+nIsh(js)),1,rout(ip),1)
  call DaXpY_(jVert,-sign,Fock(nO+1,ja+nIsh(js)),1,rout(ip),1)
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Preciba
