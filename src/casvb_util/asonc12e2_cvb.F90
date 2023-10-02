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

subroutine asonc12e2_cvb(c,axc,sxc,nvec,nprm,civb,civbh,civbs,orbs,gjorb,gjorb2,gjorb3,cvbdet,cvb)

use casvb_global, only: ipp12e, iter12e
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
integer(kind=iwp) :: nvec, nprm
real(kind=wp) :: c(nprm,nvec), axc(nprm,nvec), sxc(nprm,nvec), civb(ndet), civbh(ndet), civbs(ndet), orbs(norb,norb), gjorb(*), &
                 gjorb2(*), gjorb3(*), cvbdet(ndetvb), cvb(nvb)
integer(kind=iwp) :: ic1, ivec
real(kind=wp), allocatable :: vec_all(:)
real(kind=wp), external :: ddot_, tim_cvb

iter12e = iter12e+1
if (ipp12e >= 2) then
  write(u6,'(/,a,i5,a,f10.3,a)') ' Davidson iteration',iter12e,' at',tim_cvb(cpu0),' CPU seconds'
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
  if (.not. strucopt) call daxpy_(nvb,c(1,ivec),cvb,1,vec_all(nprorb+1),1)
  ! (CIVB set in O12EA :)
  call cizero_cvb(civbs)
  call oneexc_cvb(civb,civbs,vec_all,.false.,0)
  call str2vbc_cvb(vec_all(nprorb+1),cvbdet)
  call vb2ciaf_cvb(cvbdet,civbs)
  call cicopy_cvb(civbs,civbh)
  call makecivbhs_cvb(civbh,civbs,orbs,gjorb,gjorb2,gjorb3)

  call ci2vbg_cvb(civbh,cvbdet)
  call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
  call fzero(vec_all,nprorb)
  call onedens_cvb(civb,civbh,vec_all,.false.,0)
  call all2free_cvb(vec_all,axc(ic1,ivec),1)
  if (.not. strucopt) axc(1,ivec) = ddot_(nvb,cvb,1,vec_all(nprorb+1),1)

  call ci2vbg_cvb(civbs,cvbdet)
  call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
  call fzero(vec_all,nprorb)
  call onedens_cvb(civb,civbs,vec_all,.false.,0)
  call all2free_cvb(vec_all,sxc(ic1,ivec),1)
  if (.not. strucopt) sxc(1,ivec) = ddot_(nvb,cvb,1,vec_all(nprorb+1),1)
end do
call mma_deallocate(vec_all)

return

end subroutine asonc12e2_cvb
