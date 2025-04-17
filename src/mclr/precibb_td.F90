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

subroutine Precibb_td(ib,is,js,nd,rout,nba,Temp1,Scr,Temp2,fockii,fockai,focki,focka,Sgn)
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
use input_mclr, only: nAsh, nBas, nIsh
use Constants, only: Four, Twelve
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ib, is, js, nd, nba
real(kind=wp), intent(inout) :: rout(*)
real(kind=wp), intent(out) :: Temp1(nBa,nBa), Temp2(nBa,nBa)
real(kind=wp), intent(_OUT_) :: Scr(*)
real(kind=wp), intent(in) :: fockii, fockai, Focki(nBa,nBa), Focka(nBa,nBa), Sgn
integer(kind=iwp) :: i, ip, jVert, kB
real(kind=wp) :: ra

!                                                                      *
!***********************************************************************
!                                                                      *
jVert = nBas(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

ip = nTri_Elem(nd)-nTri_Elem(jVert)+1
ra = Four*Sgn*(Fockii+Fockai)
call COUL(jS,jS,iS,iS,iB,iB,Temp2,Scr)
Temp1(:,:) = -Sgn*Four*Temp2(:,:)
call EXCH(js,is,js,is,ib,ib,Temp2,Scr)
Temp1(:,:) = Temp1(:,:)+Sgn*Twelve*Temp2(:,:)
i = ip
do kB=nIsh(jS)+nAsh(jS)+1,nBas(jS)
  rOut(i) = rout(i)-ra
  rOut(i:i+nBas(jS)-kB) = rOut(i:i+nBas(jS)-kB)+Temp1(kB,kB:nBas(jS))+Sgn*Four*(Focki(kB,kB:nBas(jS))+Focka(kB,kB:nBas(jS)))
  i = i+nBas(jS)-kB+1
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precibb_td
