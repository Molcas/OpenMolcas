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

subroutine symtrizcvb2_cvb(vecstr,izeta,ipermzeta)

use casvb_global, only: ndetvb, norb, nsyme, nvb, nzeta
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: vecstr(nvb)
integer(kind=iwp), intent(in) :: izeta(nsyme), ipermzeta(norb,nzeta)
integer(kind=iwp) :: isyme, izeta1
real(kind=wp), allocatable :: dvbdet(:), vecstr2(:)

call mma_allocate(dvbdet,ndetvb,label='dvbdet')
call mma_allocate(vecstr2,nvb,label='vecstr2')
izeta1 = 0
do isyme=1,nsyme
  if (izeta(isyme) /= 0) then
    izeta1 = izeta1+1
    call str2vbc_cvb(vecstr,dvbdet)
    call permvb_cvb(dvbdet,ipermzeta(:,izeta1))
    call vb2strc_cvb(dvbdet,vecstr2)
    vecstr(:) = vecstr(:)+real(izeta(isyme),kind=wp)*vecstr2(:)
  end if
end do
call mma_deallocate(dvbdet)
call mma_deallocate(vecstr2)
if (izeta1 > 0) vecstr(:) = vecstr(:)/real(2**izeta1,kind=wp)

return

end subroutine symtrizcvb2_cvb
