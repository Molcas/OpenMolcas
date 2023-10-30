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

subroutine wris_cvb(ivec,n,file_id,ioffset)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp), intent(in) :: n, ivec(n)
real(kind=wp), intent(in) :: file_id
integer(kind=iwp), intent(inout) :: ioffset
integer(kind=iwp) :: ibuf(8), ilen, ioff, nreals, nrem

nreals = n/RtoI
nrem = n-nreals*RtoI
call wris_cvb_internal(ivec,ibuf)
if (nrem == 0) then
  ioffset = ioffset+nreals
else
  ioffset = ioffset+nreals+1
end if

return

! This is to allow type punning without an explicit interface
contains

#include "intent.fh"

subroutine wris_cvb_internal(ivec,ibuf)
  integer(kind=iwp), target, intent(in) :: ivec(*)
  integer(kind=iwp), target, intent(_OUT_) :: ibuf(*)
  real(kind=wp), pointer :: buf(:), vec(:)
  if (nreals > 0) then
    call c_f_pointer(c_loc(ivec(1)),vec,[nreals])
    call wrlow_cvb(vec,nreals,file_id,ioffset)
    nullify(vec)
  end if
  if (nrem > 0) then
    ilen = -1
    if (ilen >= 1+nreals+ioffset) then
      call c_f_pointer(c_loc(ibuf(1)),buf,[1])
      call rdlow_cvb(buf,1,file_id,nreals+ioffset)
      nullify(buf)
    end if
    ioff = nreals*RtoI
    ibuf(1:nrem) = ivec(ioff+1:ioff+nrem)
    call c_f_pointer(c_loc(ibuf(1)),buf,[1])
    call wrlow_cvb(buf,1,file_id,nreals+ioffset)
    nullify(buf)
  end if
end subroutine wris_cvb_internal

end subroutine wris_cvb
