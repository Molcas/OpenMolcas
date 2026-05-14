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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine indxab_cvb(indxa,indxb,nstra,nstrb,nsa,nsb)

use Symmetry_Info, only: Mul
use casvb_global, only: ityp, mxirrep, nalf, nbet, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nsa, nsb
integer(kind=iwp), intent(out) :: indxa(nsa), indxb(nsb), nstra(mxirrep), nstrb(mxirrep)
integer(kind=iwp) :: ia, ib, iisym, inda, indb, indx, irp
integer(kind=iwp), allocatable :: iocc(:)

call mma_allocate(iocc,norb+1,label='iocc')

nstra(:) = 0
nstrb(:) = 0
inda = 0
indb = 0
do iisym=1,mxirrep

  call loopstr0_cvb(iocc,indx,nalf,norb)
  do
    irp = 1
    do ia=1,nalf
      irp = Mul(irp,ityp(iocc(ia)))
    end do
    if (irp == iisym) then
      inda = inda+1
      nstra(iisym) = nstra(iisym)+1
      indxa(inda) = indx
    end if
    call loopstr_cvb(iocc,indx,nalf,norb)
    if (indx == 1) exit
  end do

  call loopstr0_cvb(iocc,indx,nbet,norb)
  do
    irp = 1
    do ib=1,nbet
      irp = Mul(irp,ityp(iocc(ib)))
    end do
    if (irp == iisym) then
      indb = indb+1
      nstrb(iisym) = nstrb(iisym)+1
      indxb(indb) = indx
    end if
    call loopstr_cvb(iocc,indx,nbet,norb)
    if (indx == 1) exit
  end do

end do

call mma_deallocate(iocc)

return

end subroutine indxab_cvb
