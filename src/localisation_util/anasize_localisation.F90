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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************

subroutine Anasize_Localisation(Den,CMO,XMO,nShell,nOrb,iSym)
! Author: T.B. Pedersen
!
! Purpose: sparsity analysis of shell-based matrices.

implicit real*8(a-h,o-z)
real*8 Den(nShell,nShell), CMO(nShell,nOrb), XMO(nShell,nOrb)
#include "WrkSpc.fh"

character*17 XHead
character*20 CHead
character*36 DHead

! Return if nothing to do.
! ------------------------

if (nShell < 0) return

! Set up bins.
! ------------

lBin = 9
call GetMem('Bin','Allo','Real',ipBin,lBin)
StpSiz = 1.0d-1
Work(ipBin) = 1.0d0
ip0 = ipBin-1
do iBin=2,lBin
  Work(ip0+iBin) = Work(ip0+iBin-1)*StpSiz
end do

! Density.
! --------

lDLT = nShell*(nShell+1)/2
call GetMem('LTDen','Allo','Real',ipDLT,lDLT)
call Sq2Tri(Den,Work(ipDLT),nShell)
write(DHead,'(A34,I2)') 'Histogram of density matrix , sym.',iSym
call Cho_Head(DHead,'=',80,6)
call Cho_Anasize(Work(ipDLT),lDlt,Work(ipBin),lBin,6)
call GetMem('LTDen','Free','Real',ipDLT,lDLT)

if (nOrb < 1) Go To 1 ! return after de-allocation

! Original MOs.
! -------------

write(CHead,'(A18,I2)') 'Original MOs, sym.',iSym
call Cho_Head(CHead,'=',80,6)
do i=1,nOrb
  write(6,'(/,2X,A,I5)') 'Original MO no.',i
  call Cho_Anasize(CMO(1,i),nShell,Work(ipBin),lBin,6)
end do

! Local MOs.
! ----------

write(XHead,'(A15,I2)') 'Local MOs, sym.',iSym
call Cho_Head(XHead,'=',80,6)
do i=1,nOrb
  write(6,'(/,2X,A,I5)') 'Local MO no.',i
  call Cho_Anasize(XMO(1,i),nShell,Work(ipBin),lBin,6)
end do

! De-allocate and return.
! -----------------------

1 call GetMem('Bin','Free','Real',ipBin,lBin)

end subroutine Anasize_Localisation
