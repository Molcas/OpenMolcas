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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nShell, nOrb, iSym
real(kind=wp), intent(in) :: Den(nShell,nShell), CMO(nShell,nOrb), XMO(nShell,nOrb)
integer(kind=iwp), parameter :: lBin = 9
integer(kind=iwp) :: i, iBin, lDLT
real(kind=wp) :: Bin(lBin), StpSiz
character(len=36) :: DHead
character(len=20) :: CHead
character(len=17) :: XHead
real(kind=wp), allocatable :: DLT(:)

! Return if nothing to do.
! ------------------------

if (nShell < 0) return

! Set up bins.
! ------------

StpSiz = 0.1_wp
Bin(1) = One
do iBin=2,lBin
  Bin(iBin) = Bin(iBin-1)*StpSiz
end do

! Density.
! --------

lDLT = nShell*(nShell+1)/2
call mma_allocate(DLT,lDLT,label='LTDen')
call Sq2Tri(Den,DLT,nShell)
write(DHead,'(A34,I2)') 'Histogram of density matrix , sym.',iSym
call Cho_Head(DHead,'=',80,u6)
call Cho_Anasize(DLT,lDlt,Bin,lBin,u6)
call mma_deallocate(DLT)

if (nOrb >= 1) then

  ! Original MOs.
  ! -------------

  write(CHead,'(A18,I2)') 'Original MOs, sym.',iSym
  call Cho_Head(CHead,'=',80,u6)
  do i=1,nOrb
    write(u6,'(/,2X,A,I5)') 'Original MO no.',i
    call Cho_Anasize(CMO(1,i),nShell,Bin,lBin,u6)
  end do

  ! Local MOs.
  ! ----------

  write(XHead,'(A15,I2)') 'Local MOs, sym.',iSym
  call Cho_Head(XHead,'=',80,u6)
  do i=1,nOrb
    write(u6,'(/,2X,A,I5)') 'Local MO no.',i
    call Cho_Anasize(XMO(1,i),nShell,Bin,lBin,u6)
  end do

end if

end subroutine Anasize_Localisation
