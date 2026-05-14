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

subroutine cnfsort_cvb(iconfs,nconf1,nel1)

use casvb_global, only: noe, norb
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nconf1, nel1
integer(kind=iwp), intent(inout) :: iconfs(noe,nconf1)
integer(kind=iwp) :: iconf, ion, iorb, jconf, mnion1, mxion1
integer(kind=iwp), allocatable :: iconfs2(:,:), ioncty(:)

call mma_allocate(ioncty,nconf1,label='ioncty')
call mma_allocate(iconfs2,noe,nconf1,label='iconfs2')

mnion1 = nel1/2
mxion1 = 0
do iconf=1,nconf1
  ion = 0
  do iorb=1,norb
    if (iconfs(iorb,iconf) == 2) ion = ion+1
  end do
  ioncty(iconf) = ion
  if (ion < mnion1) mnion1 = ion
  if (ion > mxion1) mxion1 = ion
end do
jconf = 0
do ion=mnion1,mxion1
  do iconf=1,nconf1
    if (ioncty(iconf) == ion) then
      jconf = jconf+1
      iconfs2(:,jconf) = iconfs(:,iconf)
    end if
  end do
end do
if (jconf /= nconf1) then
  write(u6,*) ' Error in cnfsort - jconf not same as nconf1 :',jconf,nconf1
  call abend_cvb()
end if
iconfs(:,:) = iconfs2(:,:)

call mma_deallocate(ioncty)
call mma_deallocate(iconfs2)

return

end subroutine cnfsort_cvb
