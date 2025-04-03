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

use Index_Functions, only: iTri, nTri_Elem
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

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = nTri_Elem(nd)
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
  ip = nTri-iTri(nd-ja+1,jVert)
  do jB=1,nAsh(jS)
    jBB = jB+nIsh(jS)
    ! Get D_(ja,jb)
    rDens = -sign*G1t((iTri(jA+nA(jS),jB+nA(jS))))
    if (jA == jB) rDens = rdens+sign*Two

    ivB = (jBB-1)*nBas(jS)+nO
    rOut(ip+1:ip+jVert) = rOut(ip+1:ip+jVert)+rDens*(Six*A_K(ivB+1:ivB+jVert)-Two*A_J(ivB+1:ivB+jVert))
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
do jA=1,nAsh(js)
  ip = nTri-iTri(nd-ja+1,nd-nAsh(js))
  rout(ip+1:ip+jVert) = rout(ip+1:ip+jVert)+sign*(Four*(FockI(nO+1:nO+jVert,ja+nIsh(js))+FockA(nO+1:nO+jVert,ja+nIsh(js)))- &
                                                  Fock(nO+1:nO+jVert,ja+nIsh(js)))
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Preciba
