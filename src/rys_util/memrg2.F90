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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine MemRg2(iAnga,nRys,MemPrm,ider)
!***********************************************************************
! This routine will compute the memory requirement of Rysg2            *
! Memory requirement is per primitive!                                 *
!                                                                      *
! Modified to gradients 1991 R. Lindh, Dept. of Theoretical Chemistry, *
! University of Lund.                                                  *
! Modified to hessians 1995 Anders Bernhardsson Dept. of Theoretical   *
! Chemistry, University of Lund.                                       *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iAnga(4), ider
integer(kind=iwp), intent(out) :: nRys, MemPrm
integer(kind=iwp) :: la, lb, lB00, lB01, lB10, lc, ld, MemPrm2, MemPrm3, n2D0, n2D1, n2D2, nabcd, nabMax, ncdMax, nPAO

!iRout = 13
!iPrint = nPrint(iRout)
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
nRys = (la+lb+lc+ld+2+ider)/2

MemPrm = 0
nPAO = nTri_Elem1(la)*nTri_Elem1(lb)*nTri_Elem1(lc)*nTri_Elem1(ld)
! 1st order gradient of [ab|cd]
MemPrm = MemPrm+nPAO*9

nabMax = la+lb+2
ncdMax = lc+ld+2
nabcd = (nabMax+1)*(ncdMax+1)
!lB10 = max(min(nabMax-1,1),0)
!lB01 = max(min(ncdMax-1,1),0)
!lB00 = max(min(Min(nabMax,ncdMax),1),0)
lB10 = 1
lB01 = 1
lB00 = 1

! 2D-Integrals

n2D0 = max(nabcd,(la+3)*(lb+3)*(ncdMax+1),(la+3)*(lb+3)*(lc+3)*(ld+3))
MemPrm = MemPrm+3*nRys*n2D0

! 1st order gradient of the 2D-integrals

n2D1 = max(nabcd,(la+3)*(lb+3)*(ncdMax+1),(la+1)*(lb+1)*(lc+1)*(ld+1)*3)
n2D2 = (la+1)*(lb+1)*(lc+1)*(ld+1)*6
MemPrm = MemPrm+3*nRys*n2D1
MemPrm3 = 3*nRys*n2D2+3*nRys

! Coefficients for recurrence relations

MemPrm2 = 3*nRys+3*nRys+3*nRys*(lB10+lB01+lB00)
! Roots
MemPrm2 = MemPrm2+nRys
! The inverse of the arguments
MemPrm2 = MemPrm2+1
! Arguments
MemPrm2 = MemPrm2+1
! Expanded versions of Zeta, ZetInv, Eta, EtaInv, P and Q
MemPrm = MemPrm+10
!if (iPrint >=99) then
!  write(u6,*) ' [ab|cd] 1st grad.   :',nPAO*9
!  write(u6,*) ' 2D-integrals        :',n2D0*3*nRys
!  write(u6,*) ' 2D-integrals (1st)  :',n2D1*3*nRys
!  write(u6,*) ' 2D-integrals (2nd)  :',n2D2*3*nRys
!  write(u6,*) ' PAQP vector         :',3*nRys
!  write(u6,*) ' QCPQ vector         :',3*nRys
!  write(u6,*) ' B10 coefficients    :',nRys*3*lB10
!  write(u6,*) ' B00 coefficients    :',nRys*3*lB00
!  write(u6,*) ' B01 coefficients    :',nRys*3*lB01
!  write(u6,*) ' Roots               :',nRys
!  write(u6,*) ' Inverse arguments   :',1
!  write(u6,*) ' Arguments           :',1
!end if
MemPrm = MemPrm+max(MemPrm2,MemPrm3)

return

end subroutine MemRg2
