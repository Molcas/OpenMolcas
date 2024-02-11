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

subroutine make_close_cvb(it)

use wadr, only: TUVX, FockOcc, DSPN, DMAT, PMAT, PA, FI, FA, D1I, D1A, OccN, CMO, DIAF, FMO
use casvb_global, only: variat
use Definitions, only: iwp
use stdalloc, only: mma_deallocate

implicit none
integer(kind=iwp) :: it
integer(kind=iwp) :: i, il, n
character(len=8) :: vec(11)
integer(kind=iwp), external :: find_lu

vec(1) = 'TMP01'
vec(2) = 'TMP02'
vec(3) = 'TMP03'
vec(4) = 'TMP04'
vec(5) = 'TMP05'
vec(6) = 'TMP06'
vec(7) = 'TMP07'
vec(8) = 'TMP08'
vec(9) = 'TMP09'
vec(10) = 'VBWFN'
vec(11) = 'JOBIPH'
il = 10
if (it == 0) il = 10
if (it == 1) il = 11
! Preassign some file names to identifiers :
do i=1,il
  n = find_lu(vec(i))
  if (n > 0) call daclos(n)
end do
if (.not. variat) then
  call mkguga_free()
  call mma_deallocate(FMO)
  call mma_deallocate(TUVX)
  Call mma_deallocate(DMAT)
  Call mma_deallocate(DSPN)
  Call mma_deallocate(PMAT)
  Call mma_deallocate(PA)
  Call mma_deallocate(DIAF)
  call mma_deallocate(FockOcc)
  call mma_deallocate(FI)
  call mma_deallocate(FA)
  call mma_deallocate(D1I)
  call mma_deallocate(D1A)
  call mma_deallocate(OccN)
  call mma_deallocate(CMO)
end if

end subroutine make_close_cvb
