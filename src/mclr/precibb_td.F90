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

subroutine Precibb_td(ib,is,js,nd,rout,nba,Temp1,Scr,Temp2,fockii,fockai,focki,focka,sign)
!***********************************************************************
!                                       [2]
! Calculates the diagonal submatrix of E    that couple
!
! kappa           with   kappa                for a
!      kinactive,virtual        kinactive,virtual
!
! single inactive index.
! Used for preconditioner.
!
! See Olsen,Yeager, Joergensen:
!  "Optimization and characterization of an MCSCF state"
!
! Called by prec
!
! ib,is       :       inactive index for the submatrix
! js          :       symmetry of virtual,virtual
! rOut        :       Submatrix
!
!***********************************************************************

use Index_Functions, only: nTri_Elem
use input_mclr, only: nAsh, nIsh, nBas
use Constants, only: Four, Twelve

implicit none
integer ib, is, js, nd
real*8 rout(*)
integer nba
real*8 Temp1(nBa,nBa)
real*8 Temp2(nBa,nBa), Scr(*)
real*8 fockii, fockai
real*8 Focki(nBa,nBa), Focka(nBa,nBa)
real*8 sign
integer jVert, ip, kB, lB
real*8 ra
integer i

!                                                                      *
!***********************************************************************
!                                                                      *
jVert = nBas(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

ip = nTri_Elem(nd)-nTri_Elem(jVert)+1
ra = Four*sign*(Fockii+Fockai)
call COUL(jS,jS,iS,iS,iB,iB,Temp2,Scr)
Temp1(:,:) = -sign*Four*Temp2(:,:)
call EXCH(js,is,js,is,ib,ib,Temp2,Scr)
Temp1(:,:) = Temp1(:,:)+sign*Twelve*Temp2(:,:)
i = ip-1
do kB=nIsh(jS)+nAsh(jS)+1,nBas(jS)
  rOut(i+1) = rout(i+1)-ra
  do lB=kb,nBAS(JS)
    i = i+1
    rOut(i) = rout(i)+Temp1(kb,lb)+sign*Four*Focki(kb,lb)+sign*Four*Focka(kb,lb)
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precibb_td
