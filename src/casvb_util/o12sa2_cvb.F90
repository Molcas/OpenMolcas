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

subroutine o12sa2_cvb(nprm,civb,civbs,cvbdet,cvb)

use casvb_global, only: ndet, ndetvb, npr, nprorb, nvb, strucopt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nprm
real(kind=wp), intent(inout) :: civb(0:ndet), civbs(0:ndet)
real(kind=wp), intent(out) :: cvbdet(ndetvb)
real(kind=wp), intent(_IN_) :: cvb(nvb)
integer(kind=iwp) :: ic1
real(kind=wp) :: cnrm, dum(1), ret
real(kind=wp), allocatable :: c(:), sxc(:), vec_all(:)
real(kind=wp), external :: ddot_

! If no optimization of structure coefficients we are doing "Augmented" calc:
if (strucopt) then
  ic1 = 1
else
  ic1 = 2
end if

call str2vbc_cvb(cvb,cvbdet)
call vb2cic_cvb(cvbdet,civb)

call cidot_cvb(civb,civbs,ret)
call ci2vbg_cvb(civbs,cvbdet)
call mma_allocate(vec_all,npr,label='vec_all')
call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
vec_all(1:nprorb) = Zero
call onedens_cvb(civb,civbs,vec_all,.false.,0)
call mma_allocate(sxc,nprm,label='sxc')
call all2free_cvb(vec_all,sxc(ic1),1)
if (.not. strucopt) sxc(1) = ddot_(nvb,cvb,1,vec_all(nprorb+1),1)
vec_all(1:nprorb) = Zero
vec_all(nprorb+1:nprorb+nvb) = cvb(:)
call mma_allocate(c,nprm,label='c')
call all2free_cvb(vec_all,c(ic1),1)
if (.not. strucopt) c(1) = ddot_(nvb,cvb,1,vec_all(nprorb+1:),1)
call mma_deallocate(vec_all)
cnrm = ddot_(nprm,c,1,sxc,1)
c(:) = c(:)/sqrt(cnrm)
sxc(:) = sxc(:)/sqrt(cnrm)
call ddrestv_cvb(c,dum,sxc,nprm,0,.false.,.true.)
call mma_deallocate(sxc)
call mma_deallocate(c)

return

end subroutine o12sa2_cvb
