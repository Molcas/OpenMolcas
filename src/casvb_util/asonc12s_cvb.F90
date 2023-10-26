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

subroutine asonc12s_cvb( &
#                       define _CALLING_
#                       include "ddasonc_interface.fh"
                       )
! Applies S and H on c vector(s).

use casvb_global, only: civb2, civb3, cpu0, cvb, cvbdet, ipp12s, iter12s, npr, nprorb, nvb, orbs, strucopt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "ddasonc_interface.fh"
integer(kind=iwp) :: ic1, ivec
real(kind=wp), allocatable :: vec_all(:)
real(kind=wp), external :: ddot_, tim_cvb

#include "macros.fh"
unused_var(axc)

iter12s = iter12s+1
if (ipp12s >= 2) then
  write(u6,'(/,a,i5,a,f10.3,a)') ' Davidson iteration',iter12s,' at',tim_cvb(cpu0),' CPU seconds'
  write(u6,'(a)') ' -----------------------------------------------'
end if

! If no optimization of structure coefficients we are doing "Augmented" calc:
if (strucopt) then
  ic1 = 1
else
  ic1 = 2
end if

call mma_allocate(vec_all,npr,label='vec_all')
do ivec=1,nvec
  call free2all_cvb(c(ic1,ivec),vec_all,1)
  if (.not. strucopt) vec_all(nprorb+1:nprorb+nvb) = vec_all(nprorb+1:nprorb+nvb)+c(1,ivec)*cvb(1:nvb)
  ! (CIVB set in O12SA :)
  call cizero_cvb(civb2)
  call oneexc_cvb(civb3,civb2,vec_all,.false.,0)
  call str2vbc_cvb(vec_all(nprorb+1),cvbdet)
  call vb2ciaf_cvb(cvbdet,civb2)
  call applyts_cvb(civb2,orbs)

  call ci2vbg_cvb(civb2,cvbdet)
  call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
  vec_all(1:nprorb) = Zero
  call onedens_cvb(civb3,civb2,vec_all,.false.,0)
  call all2free_cvb(vec_all,sxc(ic1,ivec),1)
  if (.not. strucopt) sxc(1,ivec) = ddot_(nvb,cvb,1,vec_all(nprorb+1:nprorb+nvb),1)
end do
call mma_deallocate(vec_all)

return

end subroutine asonc12s_cvb
