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

subroutine wri_cvb(ivec,n,file_id,ioffset)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: n, ivec(n), ioffset
real(kind=wp) :: file_id
#include "idbl_cvb.fh"
integer(kind=iwp) :: ibuf(8) = 0, ilen, nreals, nrem

nreals = n/idbl
nrem = n-nreals*idbl
call wri_cvb_internal(ivec,ibuf)

return

! This is to allow type punning without an explicit interface
contains

subroutine wri_cvb_internal(ivec,ibuf)
  integer(kind=iwp), target :: ivec(*), ibuf(*)
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
    call imove_cvb(ivec(1+nreals*idbl),ibuf,nrem)
    ! Trying for a "clean" write (unwritten integers written as zeros),
    ! but explicit zeroing of ibuf is not necessary as long as idbl <= 2.
    !call izero(ibuf(nrem+1),idbl-nrem)
    call c_f_pointer(c_loc(ibuf(1)),buf,[1])
    call wrlow_cvb(buf,1,file_id,nreals+ioffset)
    nullify(buf)
  end if
end subroutine wri_cvb_internal

end subroutine wri_cvb
