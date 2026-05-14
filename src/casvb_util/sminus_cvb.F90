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

subroutine sminus_cvb(bikfrom,bikto,nel,nalffrom,nalfto,nvec)

use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
real(kind=wp), intent(inout) :: bikfrom(*)
real(kind=wp), intent(_OUT_) :: bikto(*)
integer(kind=iwp), intent(in) :: nel, nalffrom, nalfto, nvec
integer(kind=iwp) :: ialffrom, ialfto, ivec, ndetfrom, ndetto
real(kind=wp) :: cnrmfrom, cnrmto
real(kind=wp), allocatable :: b4(:,:), b5(:,:)
real(kind=wp), external :: dnrm2_

call asc2ab_cvb(bikfrom,nvec,nel,nalffrom)

do ialfto=nalffrom-1,nalfto,-1
  ialffrom = ialfto+1
  call icomb_cvb(nel,ialffrom,ndetfrom)
  call icomb_cvb(nel,ialfto,ndetto)
  if (nalffrom == nalfto+1) then
    call sminus2_cvb(bikfrom,bikto,nel,ialffrom,ndetfrom,ialfto,ndetto,nvec)
  else if (ialfto == nalffrom-1) then
    call mma_allocate(b4,ndetto,nvec,label='b4')
    call sminus2_cvb(bikfrom,b4,nel,ialffrom,ndetfrom,ialfto,ndetto,nvec)
  else if (ialfto == nalfto) then
    call sminus2_cvb(b4,bikto,nel,ialffrom,ndetfrom,ialfto,ndetto,nvec)
    call mma_deallocate(b4)
  else
    call mma_allocate(b5,ndetto,nvec,label='b5')
    call sminus2_cvb(b4,b5,nel,ialffrom,ndetfrom,ialfto,ndetto,nvec)
    call mma_deallocate(b4)
    call move_alloc(b5,b4)
  end if
end do

call asc2ab_cvb(bikto,nvec,nel,nalfto)
! Now try to retain normalization ...
call icomb_cvb(nel,nalffrom,ndetfrom)
call icomb_cvb(nel,nalfto,ndetto)
do ivec=1,nvec
  cnrmfrom = dnrm2_(ndetfrom,bikfrom((ivec-1)*ndetfrom+1:ivec*ndetfrom),1)
  cnrmto = dnrm2_(ndetto,bikto((ivec-1)*ndetto+1:ivec*ndetto),1)
  if (cnrmto > 1.0e-10_wp) bikto((ivec-1)*ndetto+1:ivec*ndetto) = cnrmfrom/cnrmto*bikto((ivec-1)*ndetto+1:ivec*ndetto)
end do

return

end subroutine sminus_cvb
