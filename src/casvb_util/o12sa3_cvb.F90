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

subroutine o12sa3_cvb(vec,cvb,orbs,civec,civecp,civb,cvbdet,nparm1)

use casvb_global, only: ndet, ndetvb, norb, nprorb, nvb, strucopt
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: nparm1
real(kind=wp), intent(out) :: vec(nparm1), cvbdet(ndetvb)
real(kind=wp), intent(_IN_) :: cvb(nvb)
real(kind=wp), intent(in) :: orbs(norb,norb)
real(kind=wp), intent(inout) :: civec(0:ndet), civecp(0:ndet), civb(0:ndet)
integer(kind=iwp) :: ic1
real(kind=wp), allocatable :: vec_all(:)
real(kind=wp), external :: ddot_

call makegjorbs_cvb(orbs)

call str2vbc_cvb(cvb,cvbdet)
call vb2cic_cvb(cvbdet,civb)

call makecivecp_cvb(civec,civecp,orbs)
call ci2vbg_cvb(civecp,cvbdet)
call mma_allocate(vec_all,nparm1,label='vec_all')
call vb2strg_cvb(cvbdet,vec_all(nprorb+1))
vec_all(1:nprorb) = Zero
call onedens_cvb(civb,civecp,vec_all,.false.,0)
! If no optimization of structure coefficients we are doing "Augmented" calc:
if (strucopt) then
  ic1 = 1
else
  ic1 = 2
end if
call all2free_cvb(vec_all,vec(ic1),1)
if (.not. strucopt) vec(1) = ddot_(nvb,cvb,1,vec_all(nprorb+1:),1)
call mma_deallocate(vec_all)
call ddrhs_cvb(vec,nparm1,0)

call str2vbc_cvb(cvb,cvbdet)
call vb2cic_cvb(cvbdet,civb)

return

end subroutine o12sa3_cvb
