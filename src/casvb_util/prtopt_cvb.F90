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

subroutine prtopt_cvb()

use Index_Functions, only: nTri_Elem
use casvb_global, only: ioptim, istackrep, nfxvb, noptim, norb, nort, nvb, nzrvb, recinp
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: ioffs, ioptstep1, ioptstep2, italter, kk2, mxalter, nc_zeroed, nconvinone
integer(kind=iwp), allocatable :: idelstr(:), ifxorb(:), ifxstr(:), iorts(:,:)
logical(kind=iwp), external :: istkprobe_cvb

! First determine if end of multi-step optimization may have been reached:
if (istkprobe_cvb(istackrep)) then
  call istkpop_cvb(istackrep,nc_zeroed)
  call istkpop_cvb(istackrep,nconvinone)
  call istkpop_cvb(istackrep,italter)
  call istkpop_cvb(istackrep,mxalter)
  call istkpop_cvb(istackrep,kk2)
  call istkpop_cvb(istackrep,ioptstep2)
  call istkpop_cvb(istackrep,ioptstep1)
  call istkpush_cvb(istackrep,ioptstep1)
  call istkpush_cvb(istackrep,ioptstep2)
  call istkpush_cvb(istackrep,kk2)
  call istkpush_cvb(istackrep,mxalter)
  call istkpush_cvb(istackrep,italter)
  call istkpush_cvb(istackrep,nconvinone)
  call istkpush_cvb(istackrep,nc_zeroed)
else
  ioptstep1 = 0
  italter = 0
end if

call mma_allocate(iorts,2,nTri_Elem(norb-1),label='iorts')
call mma_allocate(ifxorb,norb,label='ifxorb')
call mma_allocate(ifxstr,nvb,label='ifxstr')
call mma_allocate(idelstr,nvb,label='idelstr')

call rdioff_cvb(11,recinp,ioffs)
call rdis_cvb(ifxorb,norb,recinp,ioffs)
call rdis_cvb(ifxstr,nfxvb,recinp,ioffs)
call rdis_cvb(idelstr,nzrvb,recinp,ioffs)
call rdis_cvb(iorts,2*nort,recinp,ioffs)

call prtopt2_cvb(ioptstep1,ioptim,italter,noptim,iorts,ifxorb,ifxstr,idelstr)

call mma_deallocate(iorts)
call mma_deallocate(ifxorb)
call mma_deallocate(ifxstr)
call mma_deallocate(idelstr)

return

end subroutine prtopt_cvb
