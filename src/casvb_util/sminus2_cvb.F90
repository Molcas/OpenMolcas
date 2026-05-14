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

subroutine sminus2_cvb(bikfrom,bikto,nel,nalffrom,ndetfrom,nalfto,ndetto,nvec)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nel, nalffrom, ndetfrom, nalfto, ndetto, nvec
real(kind=wp), intent(in) :: bikfrom(ndetfrom,nvec)
real(kind=wp), intent(out) :: bikto(ndetto,nvec)
integer(kind=iwp) :: iexc, indfrom, indto
integer(kind=iwp), allocatable :: ioccfrom(:), ioccto(:), xdetto(:,:)
integer(kind=iwp), external :: minind_cvb

call mma_allocate(xdetto,[0,nel],[0,nalfto],label='xdetto')

bikto(:,:) = Zero

! Determinant (to) weight array:
call weightfl_cvb(xdetto,nalfto,nel)
if (ndetto /= xdetto(nel,nalfto)) then
  write(u6,*) ' Discrepancy in NDET:',ndetto,xdetto(nel,nalfto)
  call abend_cvb()
end if

call mma_allocate(ioccfrom,nalffrom,label='ioccfrom')
call mma_allocate(ioccto,nalfto,label='ioccto')

call loopstr0_cvb(ioccfrom,indfrom,nalffrom,nel)
do
  ioccto(1:nalfto) = ioccfrom(2:nalfto+1)
  do iexc=1,nalffrom
    indto = minind_cvb(ioccto,nalfto,nel,xdetto)
    bikto(indto,:) = bikto(indto,:)+bikfrom(indfrom,:)
    if (iexc < nalffrom) ioccto(iexc) = ioccfrom(iexc)
  end do
  call loopstr_cvb(ioccfrom,indfrom,nalffrom,nel)
  if (indfrom == 1) exit
end do

call mma_deallocate(xdetto)
call mma_deallocate(ioccfrom)
call mma_deallocate(ioccto)

return

end subroutine sminus2_cvb
