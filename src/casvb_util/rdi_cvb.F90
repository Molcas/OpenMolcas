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

subroutine rdi_cvb(ivec,n,file_id,ioffset)

use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
use Definitions, only: wp, iwp, RtoI

implicit none
integer(kind=iwp) :: n, ivec(n), ioffset
real(kind=wp) :: file_id
integer(kind=iwp) :: ibuf(8), nreals, nrem

nreals = n/RtoI
nrem = n-nreals*RtoI
call rdi_cvb_internal(ivec,ibuf)

return

! This is to allow type punning without an explicit interface
contains

subroutine rdi_cvb_internal(ivec,ibuf)
  integer(kind=iwp), target :: ivec(*), ibuf(*)
  real(kind=wp), pointer :: buf(:), vec(:)
  call c_f_pointer(c_loc(ivec(1)),vec,[nreals])
  call rdlow_cvb(vec,nreals,file_id,ioffset)
  nullify(buf)
  if (nrem > 0) then
    call c_f_pointer(c_loc(ibuf(1)),buf,[1])
    call rdlow_cvb(buf,1,file_id,nreals+ioffset)
    nullify(buf)
    call imove_cvb(ibuf,ivec(1+nreals*RtoI),nrem)
  end if
end subroutine rdi_cvb_internal

end subroutine rdi_cvb
