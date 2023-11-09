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

subroutine rdcivec_cvb(detvec,fn,reord)
!***********************************************************************
!                                                                      *
!     Read the contents of the JOBIPH file.                            *
!                                                                      *
!***********************************************************************

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: detvec(*)
character(len=*), intent(in) :: fn
logical(kind=iwp) :: reord
#include "rasdim.fh"
#include "jobiph_j.fh"
integer(kind=iwp) :: i, iDisk, iwr, j, k, lujob, ncix(8), ndet_j
real(kind=wp) :: rdum(1)
real(kind=wp), allocatable :: CI(:), CI2(:)
logical(kind=iwp), parameter :: debug = .false.

iwr = 0

call getnci_cvb(ncix,nactel_j,ispin_j-1,lsym_j)
ndet_j = ncix(1)

lujob = 15
call daname_cvb(lujob,fn)
! Allocate at least NDET words for each vector, since this is
! required by csdtvc:
!call mma_allocate(CI,nConf_j*nroots_j,label='OCIvec')
call mma_allocate(CI,nConf_j*nroots_j+ndet_j-nConf_j,label='OCIvec')
if (iwr == 0) then
  do i=1,nroots_j
    j = iroot_j(i)
    iDisk = iadr15_j(4)
    do k=1,j-1
      call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
    end do
    call dDaFile(LuJob,2,CI(1+(i-1)*nConf_j),nConf_j,iDisk)
  end do

  if (reord) then
    call mma_allocate(CI2,nConf_j,label='CI2')
    call reord2_cvb(CI,CI2,1)
    CI(1:nConf_j) = CI2(:)
    call mma_deallocate(CI2)
  end if

  call csf2det_cvb(CI,detvec,lsym_j,1)
else if (iwr == 1) then
  call csf2det_cvb(CI,detvec,lsym_j,2)

  if (reord) then
    call mma_allocate(CI2,nConf_j,label='CI2')
    call reord2_cvb(CI,CI2,0)
    CI(1:nConf_j) = CI2(:)
    call mma_deallocate(CI2)
  end if

  do i=1,nroots_j
    j = iroot_j(i)
    iDisk = iadr15_j(4)
    do k=1,j-1
      call dDaFile(LuJob,0,rdum,nConf_j,iDisk)
    end do
    call dDaFile(LuJob,1,CI(1+(i-1)*nConf_j),nConf_j,iDisk)
  end do
end if
if (debug) then
  do i=0,nroots_j-1
    write(u6,'(a,i3,a)') ' (CSF) CI vector ',i+1,' :'
    write(u6,'(a)') ' ---------------------'
    call mxprint_cvb(CI(1+nconf_j*i),1,nconf_j,0)
  end do
end if
call mma_deallocate(CI)
call daclos_cvb(lujob)

return

end subroutine rdcivec_cvb
