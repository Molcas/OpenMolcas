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

subroutine Precibb(ib,is,js,nd,rout,no,Temp1,Scr,Temp2,fockii,fockai,focki,focka,sign)
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
use input_mclr, only: nAsh, nIsh, nOrb
use Constants, only: Four, Twelve

implicit none
integer ib, is, js, nd
real*8 rout(*)
integer no
real*8 Temp2(*), Temp1(*), Scr(*)
real*8 fockii, fockai
real*8 Focki(no,no), Focka(no,no)
real*8 sign
integer jVert, ip, kB, lB
real*8 ra
integer i

!                                                                      *
!***********************************************************************
!                                                                      *
jVert = nOrb(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

ip = nTri_Elem(nd)-nTri_Elem(jVert)+1
ra = Four*sign*(Fockii+Fockai)
call COUL(jS,jS,iS,is,IB,iB,Temp2,Scr)
call Dyax(no**2,-sign*Four,Temp2,1,Temp1,1)
call EXCH(js,is,js,is,ib,ib,Temp2,Scr)
call DaXpY_(no**2,sign*Twelve,Temp2,1,Temp1,1)
i = ip-1
do kB=nIsh(jS)+nAsh(jS),nOrb(jS)-1
  rOut(i+1) = rout(i+1)-ra
  do lB=kb,nOrb(JS)-1
    i = i+1
    rOut(i) = rout(i)+Temp1(kb+1+no*lb)+sign*Four*Focki(kb+1,lb+1)+sign*Four*Focka(kb+1,lb+1)
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Precibb
