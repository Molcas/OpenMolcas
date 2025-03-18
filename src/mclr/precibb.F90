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

subroutine Precibb(ib,is,js,nd,rout,nba,no,Temp1,Scr,Temp2,fockii,fockai,focki,focka,sign)
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

use input_mclr, only: nAsh, nIsh, nOrb

implicit none
integer ib, is, js, nd
real*8 rout(*)
integer nba, no
real*8 Temp2(*), Temp1(*), Scr(*)
real*8 fockii, fockai
real*8 Focki(no,no), Focka(no,no)
real*8 sign
integer nTri, jVert, i1, ip, kB, lB
real*8 ra
! Statement functions
integer i, j, iTri, iTri1
iTri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
iTri1(i,j) = nTri-itri(nd-min(i,j)+1,nd-min(i,j)+1)+max(i,j)-min(i,j)+1

!                                                                      *
!***********************************************************************
!                                                                      *
nTri = itri(nd,nd)

jVert = nOrb(js)-nAsh(js)-nIsh(js)
if (jvert == 0) return

i1 = nD-jVert+1
ip = itri1(i1,i1)
ra = 4.0d0*sign*(Fockii+Fockai)
call COUL(jS,jS,iS,is,IB,iB,Temp2,Scr)
call Dyax(no**2,-sign*4.0d0,Temp2,1,Temp1,1)
call EXCH(js,is,js,is,ib,ib,Temp2,Scr)
call DaXpY_(no**2,sign*12.0d0,Temp2,1,Temp1,1)
i = ip-1
do kB=nIsh(jS)+nAsh(jS),nOrb(jS)-1
  rOut(i+1) = rout(i+1)-ra
  do lB=kb,nOrb(JS)-1
    i = i+1
    rOut(i) = rout(i)+Temp1(kb+1+no*lb)+sign*4.0d0*Focki(kb+1,lb+1)+sign*4.0d0*Focka(kb+1,lb+1)
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Avoid unused argument warnings
if (.false.) call Unused_integer(nba)

end subroutine Precibb
