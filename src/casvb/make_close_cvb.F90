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

use wadr, only: TUVX
use casvb_global, only: ipfocc_cvb, lcmo_cvb, ld1a_cvb, ld1i_cvb, ld1tot_cvb, ldiaf_cvb, ldmat_cvb, ldspn_cvb, lfa_cvb, lfi_cvb, &
                        loccn_cvb, lpa_cvb, lpmat_cvb, lw1_cvb, variat
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
  call getmem('CICTL1','FREE','REAL',lw1_cvb,0)
  call mma_deallocate(TUVX)
  call getmem('DMAT','FREE','REAL',ldmat_cvb,0)
  call getmem('DSPN','FREE','REAL',ldspn_cvb,0)
  call getmem('PMAT','FREE','REAL',lpmat_cvb,0)
  call getmem('P2AS','FREE','REAL',lpa_cvb,0)
  call getmem('DIAF','FREE','REAL',ldiaf_cvb,0)
  call getmem('FOCC','FREE','REAL',ipfocc_cvb,0)
  call getmem('FI','FREE','REAL',lfi_cvb,0)
  call getmem('FA','FREE','REAL',lfa_cvb,0)
  call getmem('D1I','FREE','REAL',ld1i_cvb,0)
  call getmem('D1A','FREE','REAL',ld1a_cvb,0)
  call getmem('D1tot','FREE','REAL',ld1tot_cvb,0)
  call getmem('OCCN','FREE','REAL',loccn_cvb,0)
  call getmem('LCMO','FREE','REAL',lcmo_cvb,0)
end if

return

end subroutine make_close_cvb
